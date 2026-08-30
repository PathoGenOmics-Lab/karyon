// A camera over a tree that does not move.
//
// This is the other way to look at a phylogeny, and the one every viewer of a
// million tips uses. The crate works the layout out once and hands over the
// coordinates; nothing here decides where a branch goes. What moves is the
// window onto them, so a wheel notch costs a repaint of what is on screen and
// never a walk over the tree, whatever the tree holds.
//
// Two things keep a repaint cheap however big the file is.
//
// The rows are a coordinate, so sorting the nodes by row once makes the rows in
// view a contiguous run: a binary search at each edge of the window and the
// work is proportional to what can be seen and not to what exists.
//
// And a screen has a few hundred rows of pixels whatever the tree has of rows.
// Past one row per pixel the extra ones land on top of each other, so they are
// skipped: the same reason every track in this crate bins to a point per pixel
// column, applied to the other axis.

window.karyonCanvas = (function () {
  "use strict";

  var LABEL_ROOM = 9;

  // The rail down the right, and the narrowest canvas that gets one. It is a
  // scrollbar that shows what it is scrolling through, so it is read at a
  // glance and never zoomed: it always holds the whole tree.
  var RAIL_WIDE = 78;
  var RAIL_NARROW = 34;
  var RAIL_PAD = 5;
  // Narrow enough that a phone keeps one. It used to stand down below 420,
  // which is every phone there is, so the reader with the least screen to know
  // where they were in got the least help finding out.
  var RAIL_LEAST = 260;
  // A window of three rows in two million is a thousandth of a pixel tall. The
  // mark for it stays this big so there is always something to see and to grab.
  var MARK_LEAST = 4;

  // The dial: the whole disc, small, in the corner furthest from the middle of
  // the canvas, and the gap it keeps from the edge.
  var DIAL_SIDE = 132;
  var DIAL_EDGE = 10;

  function make(canvas) {
    var view = null;
    var camera = null;
    // The camera the two flat projections share, and which of the three is
    // being looked through: "rows", "disc" or "spread".
    var disc = null;
    var mode = "rows";
    var decoder = new TextDecoder();
    var encoder = new TextEncoder();

    function nameOf(node) {
      var len = view.length[node];
      if (!len) return "";
      var at = view.start[node];
      return decoder.decode(view.names.subarray(at, at + len));
    }

    // Everything by row, once, so the rows in view are a run rather than a
    // search.
    function order(placed) {
      var by = new Uint32Array(placed.count);
      for (var i = 0; i < placed.count; i++) by[i] = i;
      var y = placed.y;
      var sorted = Array.prototype.slice.call(by);
      sorted.sort(function (a, b) { return y[a] - y[b]; });
      return Uint32Array.from(sorted);
    }

    // `keep` says this is the same tree laid out again, so the window onto it
    // still means something. Without it a reader who had navigated into a
    // clade lost it every time they pressed Cladogram.
    function load(placed, keep) {
      var was = keep && view && view.count === placed.count
        ? { camera: camera, disc: disc, mode: mode }
        : null;
      view = placed;
      // A layout with no root sends the order to read its terminals in, since
      // it has no rows for the page to sort by. That order is what stands in
      // for rows here: it is what one is picked out of every few from.
      view.found = null;
      view.tipsBeyond = null;
      view.rootless = !!(placed.order && placed.order.length);
      view.byRow = view.rootless ? placed.order : order(placed);
      mode = view.rootless ? "spread" : "rows";
      view.dial = null;
      view.overview = null;
      view.overviewFor = -1;
      var lowX = Infinity, highX = -Infinity, lowY = Infinity, highY = -Infinity;
      for (var i = 0; i < placed.count; i++) {
        if (placed.x[i] < lowX) lowX = placed.x[i];
        if (placed.x[i] > highX) highX = placed.x[i];
        if (placed.y[i] < lowY) lowY = placed.y[i];
        if (placed.y[i] > highY) highY = placed.y[i];
      }
      if (!(highX > lowX)) highX = lowX + 1;
      if (!(highY > lowY)) highY = lowY + 1;
      view.bounds = { lowX: lowX, highX: highX, lowY: lowY, highY: highY };
      // Room to mark who is already on the list while a paint walks upward,
      // and a counter so the marks are told apart by age rather than cleared.
      // The root, which every walk ends at and which the pixel test needs the
      // depth of before any walk has run.
      view.root = 0;
      for (var look = 0; look < placed.count; look++) {
        if (placed.parent[look] === 0xffffffff) { view.root = look; break; }
      }
      view.seen = new Int32Array(placed.count);
      view.visit = 0;
      view.shown = { child: new Uint32Array(1 << 14), up: new Uint32Array(1 << 14) };
      if (view.rootless) measureSpread();
      home();
      if (was && was.mode === mode) {
        // The rows are in the same order whichever way the branches are
        // measured, so the window onto them still points at the same tips. The
        // depths have moved, so the camera in x is left to be fitted again.
        if (mode === "rows") {
          camera.y0 = was.camera.y0;
          camera.y1 = was.camera.y1;
        } else if (was.disc) {
          disc.cx = was.disc.cx;
          disc.cy = was.disc.cy;
          disc.half = was.disc.half;
          settleDisc();
        }
      }
    }

    function home() {
      if (mode !== "rows") { discHome(); return; }
      var b = view.bounds;
      // A margin on the right for the names, which are drawn outward from the
      // tip and are not in the coordinates.
      camera = { x0: b.lowX, x1: b.highX + (b.highX - b.lowX) * 0.25, y0: b.lowY - 0.5, y1: b.highY + 0.5 };
    }

    function size() {
      var ratio = window.devicePixelRatio || 1;
      // A canvas with no size in the stylesheet takes its width from its own
      // backing store, so writing the store below moves the box that was just
      // measured and the two chase each other: a page that forgot the CSS grew
      // one to sixteen million pixels a side in a handful of frames. The stop
      // costs nothing on a page that sized its canvas, and turns a hang into a
      // picture that is merely wrong on one that did not.
      var most = 8192 / ratio;
      var wide = Math.min(canvas.clientWidth || 1, most);
      var tall = Math.min(canvas.clientHeight || 1, most);
      if (canvas.width !== Math.round(wide * ratio) || canvas.height !== Math.round(tall * ratio)) {
        canvas.width = Math.round(wide * ratio);
        canvas.height = Math.round(tall * ratio);
      }
      return { wide: wide, tall: tall, ratio: ratio };
    }

    // Where the rail is, or null on a canvas too narrow to give it the room.
    function rail(box) {
      if (!view || box.wide < RAIL_LEAST) return null;
      // A share of the width, held between a size worth drawing and a size
      // worth giving up for.
      var wide = Math.max(RAIL_NARROW, Math.min(RAIL_WIDE, Math.round(box.wide * 0.14)));
      return { x0: box.wide - wide, wide: wide, tall: box.tall };
    }

    // The rail holds every row there is, so a row maps to a height on it and a
    // height back to a row. This is the only arithmetic the rail needs, and it
    // is the inverse of itself.
    function rowAtHeight(py, box) {
      var b = view.bounds;
      var span = b.highY - b.lowY + 1;
      return b.lowY + (py / Math.max(1, box.tall)) * span;
    }

    function heightOfRow(row, box) {
      var b = view.bounds;
      var span = b.highY - b.lowY + 1;
      return ((row - b.lowY) / span) * box.tall;
    }

    // The whole tree, thinned to the rail's own height. Built once, and again
    // only if the canvas changes height or the page changes colour, because it
    // is the same picture at a different resolution and not a second opinion
    // about the tree: it comes out of the same `select` the main view is drawn
    // from.
    //
    // The picture never moves, so it is drawn once onto a canvas of its own and
    // stamped from there. Re-stroking it every frame was the floor under every
    // gesture: on a tree with six hundred thousand nodes it cost ten
    // milliseconds a frame even with fifty branches on screen.
    function overview(box, theme) {
      var strip = rail(box);
      var across = strip ? strip.wide : RAIL_WIDE;
      if (
        view.overviewFor === box.tall &&
        view.overviewInk === theme.faint &&
        view.overviewWide === across
      ) {
        return view.overview;
      }
      var runs = [[0, view.byRow.length]];
      var stride = Math.max(1, Math.ceil(view.byRow.length / Math.max(1, box.tall)));
      var picked = select(
        runs,
        box.tall,
        rowPixels(runs, across, box.tall, view.bounds.highY - view.bounds.lowY + 1, stride),
        { child: new Uint32Array(1 << 13), up: new Uint32Array(1 << 13) }
      );
      var seen = {
        child: picked.child,
        up: picked.up,
        count: picked.count,
        reach: spread(picked),
      };
      seen.wide = across;
      seen.plate = plate(box, theme, seen);
      view.overview = seen;
      view.overviewFor = box.tall;
      view.overviewInk = theme.faint;
      view.overviewWide = across;
      return seen;
    }

    // A canvas holding just the rail's silhouette, or null where there is no
    // document to make one from, which is where the checks run. Drawing to the
    // main canvas is the fallback and it is the same code either way.
    function plate(box, theme, seen) {
      var papers = canvas.ownerDocument;
      if (!papers || !papers.createElement) return null;
      var sheet = papers.createElement("canvas");
      var ratio = box.ratio;
      sheet.width = Math.max(1, Math.round(seen.wide * ratio));
      sheet.height = Math.max(1, Math.round(box.tall * ratio));
      var ink = sheet.getContext("2d");
      if (!ink) return null;
      ink.setTransform(ratio, 0, 0, ratio, 0, 0);
      strokeOverview(ink, theme, box, { x0: 0, wide: seen.wide }, seen);
      return sheet;
    }

    // The silhouette itself, given somewhere to put it.
    function strokeOverview(ink, theme, box, strip, seen) {
      var inner = Math.max(4, seen.wide - RAIL_PAD * 2);
      var reach = seen.reach;
      var acrossX = reach.highX - reach.lowX || 1;
      var atX = function (value) {
        return strip.x0 + RAIL_PAD + ((value - reach.lowX) / acrossX) * inner;
      };
      var atY = function (row) { return heightOfRow(row, box); };
      ink.strokeStyle = theme.faint;
      ink.lineWidth = 1;
      ink.beginPath();
      for (var each = 0; each < seen.count; each++) {
        var node = seen.child[each];
        var over = seen.up[each];
        var y = atY(view.y[node]);
        ink.moveTo(atX(view.x[over]), y);
        ink.lineTo(atX(view.x[node]), y);
        ink.moveTo(atX(view.x[over]), y);
        ink.lineTo(atX(view.x[over]), atY(view.y[over]));
      }
      ink.stroke();
    }

    // The first row at or after `value`, by binary search over the sorted rows.
    function firstRow(value) {
      var low = 0, high = view.byRow.length;
      while (low < high) {
        var mid = (low + high) >> 1;
        if (view.y[view.byRow[mid]] < value) low = mid + 1;
        else high = mid;
      }
      return low;
    }

    // The disc's own way of saying where you are. A circle has no top and
    // bottom, so the rail's bar does not fit it: what fits is the whole disc,
    // small, in a corner, with the window drawn on it. Same idea, shape
    // following the projection.
    function dialPlate(theme, side) {
      var papers = canvas.ownerDocument;
      var span = side / 2 / 1.1;
      var budget = Math.max(64, Math.round(Math.PI * 2 * span));
      var stepA = budget / (Math.PI * 2);
      // The dial reads the same pixels the view it stands for does, at its own
      // size: how far out and how far round on the disc, and plainly across and
      // down where there is no middle to measure from.
      var pixel = mode === "disc"
        ? {
            column: function (node) { return (radiusOf(view.x[node]) * span) | 0; },
            row: function (node) { return (angleOf(view.y[node]) * stepA) | 0; },
          }
        : {
            column: function (node) { return (unitX(node) * span) | 0; },
            row: function (node) { return (unitY(node) * span) | 0; },
          };
      var picked = select(
        [[0, view.byRow.length]],
        budget,
        pixel,
        { child: new Uint32Array(1 << 12), up: new Uint32Array(1 << 12) }
      );
      var held = { child: picked.child, up: picked.up, count: picked.count, span: span };
      if (!papers || !papers.createElement) return held;
      var sheet = papers.createElement("canvas");
      var ratio = window.devicePixelRatio || 1;
      sheet.width = Math.round(side * ratio);
      sheet.height = Math.round(side * ratio);
      var ink = sheet.getContext("2d");
      if (!ink) return held;
      ink.setTransform(ratio, 0, 0, ratio, 0, 0);
      strokeDisc(ink, theme, held, side / 2, side / 2, span);
      held.plate = sheet;
      return held;
    }

    function strokeDisc(ink, theme, held, midX, midY, span) {
      ink.strokeStyle = theme.faint;
      ink.lineWidth = 1;
      ink.beginPath();
      if (mode !== "disc") {
        // Straight lines between two points, because that is what the
        // projection is at any size.
        for (var line = 0; line < held.count; line++) {
          var from = held.child[line];
          var to = held.up[line];
          ink.moveTo(midX + unitX(to) * span, midY + unitY(to) * span);
          ink.lineTo(midX + unitX(from) * span, midY + unitY(from) * span);
        }
        ink.stroke();
        return;
      }
      for (var each = 0; each < held.count; each++) {
        var node = held.child[each];
        var over = held.up[each];
        var turn = angleOf(view.y[node]);
        var back = angleOf(view.y[over]);
        var inward = radiusOf(view.x[over]) * span;
        var out = radiusOf(view.x[node]) * span;
        var alongX = midX + Math.cos(turn) * inward;
        var alongY = midY + Math.sin(turn) * inward;
        if ((inward * (turn - back) * (turn - back)) / 8 < 0.5) {
          ink.moveTo(midX + Math.cos(back) * inward, midY + Math.sin(back) * inward);
          ink.lineTo(alongX, alongY);
        } else {
          ink.moveTo(midX + Math.cos(back) * inward, midY + Math.sin(back) * inward);
          ink.arc(midX, midY, inward, back, turn, turn < back);
        }
        ink.moveTo(alongX, alongY);
        ink.lineTo(midX + Math.cos(turn) * out, midY + Math.sin(turn) * out);
      }
      ink.stroke();
    }

    function dial(box) {
      if (!view || box.wide < RAIL_LEAST) return null;
      // Never more than about a third of either side, so a small canvas keeps
      // a small dial rather than losing a third of its picture to one.
      var side = Math.round(
        Math.max(64, Math.min(DIAL_SIDE, box.wide * 0.34, box.tall * 0.34))
      );
      return {
        x0: box.wide - side - DIAL_EDGE,
        y0: box.tall - side - DIAL_EDGE,
        side: side,
      };
    }

    function paintDial(ctx, theme, box, arc) {
      var spot = dial(box);
      if (!spot) return null;
      if (
        !view.dial ||
        view.dialInk !== theme.faint ||
        view.dialMode !== mode ||
        view.dialSide !== spot.side
      ) {
        view.dial = dialPlate(theme, spot.side);
        view.dialInk = theme.faint;
        view.dialMode = mode;
        view.dialSide = spot.side;
      }
      var held = view.dial;

      ctx.fillStyle = theme.plate;
      ctx.fillRect(spot.x0, spot.y0, spot.side, spot.side);
      ctx.strokeStyle = theme.frame;
      ctx.lineWidth = 1;
      ctx.strokeRect(spot.x0 + 0.5, spot.y0 + 0.5, spot.side - 1, spot.side - 1);

      var midX = spot.x0 + spot.side / 2;
      var midY = spot.y0 + spot.side / 2;
      if (held.plate) ctx.drawImage(held.plate, spot.x0, spot.y0, spot.side, spot.side);
      else strokeDisc(ctx, theme, held, midX, midY, held.span);

      // The window, on the disc. At the far end of a zoom it is a speck, so it
      // is held to a size a reader can see and a hand can catch, the way the
      // rail holds its bar.
      var seat = discBox(box);
      var wide = Math.max(MARK_LEAST, 2 * seat.halfW * held.span);
      var tall = Math.max(MARK_LEAST, 2 * disc.half * held.span);
      var left = midX + (disc.cx - seat.halfW) * held.span;
      var top = midY + (disc.cy - disc.half) * held.span;
      if (wide > spot.side) { left = spot.x0; wide = spot.side; }
      if (tall > spot.side) { top = spot.y0; tall = spot.side; }
      if (left < spot.x0) left = spot.x0;
      if (top < spot.y0) top = spot.y0;
      if (left + wide > spot.x0 + spot.side) left = spot.x0 + spot.side - wide;
      if (top + tall > spot.y0 + spot.side) top = spot.y0 + spot.side - tall;
      ctx.fillStyle = theme.window;
      ctx.fillRect(left, top, wide, tall);
      ctx.strokeStyle = theme.edge;
      ctx.strokeRect(left + 0.5, top + 0.5, Math.max(1, wide - 1), Math.max(1, tall - 1));

      return {
        x0: spot.x0, y0: spot.y0, side: spot.side,
        top: top, left: left, deep: tall, wide: wide,
        drawn: held.count,
        arc: arc ? arc.span : Math.PI * 2,
      };
    }

    // ------------------------------------------------------- the flat views

    // Both projections without rows put a node somewhere in a plane, and both
    // are looked at through the same camera: a box over a space where the
    // drawing is two units across. The disc bends the rows round; the rootless
    // walk has no rows to bend, and its coordinates arrive already in a plane.
    // Everything below this line is shared by the two of them.
    function unitX(node) {
      if (mode === "disc") return Math.cos(angleOf(view.y[node])) * radiusOf(view.x[node]);
      return (view.x[node] - view.midX) / view.reach;
    }

    function unitY(node) {
      if (mode === "disc") return Math.sin(angleOf(view.y[node])) * radiusOf(view.x[node]);
      return (view.y[node] - view.midY) / view.reach;
    }

    // What the rootless walk needs before it can be looked at: where the middle
    // of the drawing is and how far it reaches, so it lands in the same two
    // units across the disc does and the camera does not have to know which it
    // is looking at.
    function measureSpread() {
      var b = view.bounds;
      view.midX = (b.lowX + b.highX) / 2;
      view.midY = (b.lowY + b.highY) / 2;
      view.reach = Math.max(b.highX - b.lowX, b.highY - b.lowY) / 2 || 1;
    }

    // Which sides of the window a point falls outside, so an edge with both
    // ends outside the same side can be dropped without asking whether it
    // crosses. An edge with ends outside different sides may still cross the
    // window, and is kept.
    function outside(px, py, box) {
      var pad = 24;
      var code = 0;
      if (px < -pad) code |= 1;
      if (px > box.wide + pad) code |= 2;
      if (py < -pad) code |= 4;
      if (py > box.tall + pad) code |= 8;
      return code;
    }

    // ----------------------------------------------------------------- disc

    // The circular projection is the rectangular one in polar coordinates: the
    // same rows, with depth becoming radius and row becoming angle. The crate
    // says exactly that in radial.rs, and this is its arithmetic rather than
    // another one that looks similar, so the canvas and the figure agree. The
    // three constants are the crate's own defaults: the sweep starts at the top
    // of the circle, goes all the way round, and leaves a hole in the middle
    // eight percent of the way out.
    var DISC_START = -Math.PI / 2;
    var DISC_SWEEP = Math.PI * 2;
    var DISC_HOLE = 0.08;

    function terminals() {
      if (view.rootless) return Math.max(1, view.byRow.length);
      var b = view.bounds;
      return Math.max(1, b.highY - b.lowY + 1);
    }

    function angleOf(row) {
      var b = view.bounds;
      var many = terminals();
      if (many <= 1) return DISC_START;
      return DISC_START + (DISC_SWEEP * (row - b.lowY)) / many;
    }

    function rowAtAngle(angle) {
      var b = view.bounds;
      return b.lowY + ((angle - DISC_START) / DISC_SWEEP) * terminals();
    }

    // Depth as a fraction of the way from the hole to the rim, which is what
    // the crate's `radius` does with the scene's own minimum and maximum.
    function radiusOf(depth) {
      var b = view.bounds;
      var across = b.highX - b.lowX;
      var part = across > 0 ? (depth - b.lowX) / across : 0;
      if (part < 0) part = 0;
      if (part > 1) part = 1;
      return DISC_HOLE + part * (1 - DISC_HOLE);
    }

    // The disc is drawn in a space where the rim is one unit from the middle,
    // and the camera is a box over that space rather than a run of rows. A
    // circle has no top and bottom to scroll between, so what moves is a
    // window over a map.
    function discHome() {
      disc = { cx: 0, cy: 0, half: 1.1 };
    }

    function discBox(box) {
      var scale = box.tall / (2 * disc.half);
      return {
        scale: scale,
        midX: box.wide / 2,
        midY: box.tall / 2,
        halfW: disc.half * (box.wide / box.tall),
      };
    }

    // Which angles the window can see. Null means all of them, which is what a
    // window holding the middle of the disc sees however small it is.
    function arcInView(box) {
      var seat = discBox(box);
      var x0 = disc.cx - seat.halfW, x1 = disc.cx + seat.halfW;
      var y0 = disc.cy - disc.half, y1 = disc.cy + disc.half;
      if (x0 <= 0 && x1 >= 0 && y0 <= 0 && y1 >= 0) return null;
      var bearings = [];
      var steps = 24;
      for (var i = 0; i <= steps; i++) {
        var t = i / steps;
        bearings.push(Math.atan2(y0, x0 + (x1 - x0) * t));
        bearings.push(Math.atan2(y1, x0 + (x1 - x0) * t));
        bearings.push(Math.atan2(y0 + (y1 - y0) * t, x0));
        bearings.push(Math.atan2(y0 + (y1 - y0) * t, x1));
      }
      bearings.sort(function (a, b) { return a - b; });
      // The widest gap between one angle and the next is the part of the circle
      // the window cannot see, so what it can see is everything else.
      var gap = bearings[0] + Math.PI * 2 - bearings[bearings.length - 1];
      var at = bearings.length - 1;
      for (var k = 0; k + 1 < bearings.length; k++) {
        var wideEnough = bearings[k + 1] - bearings[k];
        if (wideEnough > gap) { gap = wideEnough; at = k; }
      }
      var from = bearings[(at + 1) % bearings.length];
      var to = bearings[at];
      if (to < from) to += Math.PI * 2;
      return { from: from, to: to, span: to - from };
    }

    // The rows those angles stand for, as spans of `byRow`. A window that
    // straddles the place where the circle's ends meet gets two.
    function runsForArc(arc) {
      var everything = [[0, view.byRow.length]];
      if (!arc) return everything;
      var b = view.bounds;
      var many = terminals();
      var first = rowAtAngle(arc.from) - 1;
      var last = rowAtAngle(arc.to) + 1;
      if (last - first >= many) return everything;
      if (first < b.lowY) {
        return [[firstRow(first + many), view.byRow.length], [0, firstRow(last)]];
      }
      if (last > b.highY) {
        return [[firstRow(first), view.byRow.length], [0, firstRow(last - many)]];
      }
      return [[firstRow(first), firstRow(last)]];
    }

    // The disc, drawn. A branch is the same elbow the rectangular view draws,
    // bent: an arc at the parent's radius from the parent's angle round to the
    // child's, then a straight run outward at the child's angle.
    function paintDisc(theme) {
      var box = size();
      var ctx = canvas.getContext("2d");
      ctx.setTransform(box.ratio, 0, 0, box.ratio, 0, 0);
      ctx.clearRect(0, 0, box.wide, box.tall);

      var seat = discBox(box);
      var midX = seat.midX - disc.cx * seat.scale;
      var midY = seat.midY - disc.cy * seat.scale;
      var arc = arcInView(box);
      var runs = runsForArc(arc);

      // A screen can tell apart as many angles as the widest circle it can see
      // has pixels along the part of it in view, which is the circular reading
      // of one row per pixel. The widest circle it can see is not always the
      // rim: zoomed in on the middle of the disc the rim is off the canvas
      // entirely, and measuring against it asked for six thousand samples of a
      // region that can hold a few hundred.
      var farX = Math.max(Math.abs(disc.cx - seat.halfW), Math.abs(disc.cx + seat.halfW));
      var farY = Math.max(Math.abs(disc.cy - disc.half), Math.abs(disc.cy + disc.half));
      var outermost = Math.min(1, Math.sqrt(farX * farX + farY * farY));
      var rim = (arc ? arc.span : Math.PI * 2) * outermost * seat.scale;
      var budget = Math.max(64, Math.min(4000, Math.round(rim)));
      var limits = {
        stop: inside > 0 ? function (node) { return radiusOf(view.x[node]) < inside; } : null,
        skip: outermost < 1
          ? function (walk, over) { return radiusOf(view.x[over]) > outermost; }
          : null,
      };
      // Every walk heads inward, so what ends it is being further in than
      // anything the window can see. That is the distance from the middle of
      // the disc to the nearest corner or edge of the window, and it is zero
      // when the window holds the middle. Stopping on "off the canvas" instead
      // was wrong in a way worth writing down: a tip outside the window walking
      // inward would be cut off before it reached the part of its own lineage
      // that is inside it, and the picture came out with 1,423 of its 1,430
      // branches off screen.
      var seatX0 = disc.cx - seat.halfW, seatX1 = disc.cx + seat.halfW;
      var seatY0 = disc.cy - disc.half, seatY1 = disc.cy + disc.half;
      var nearX = Math.max(seatX0, Math.min(0, seatX1));
      var nearY = Math.max(seatY0, Math.min(0, seatY1));
      var inside = Math.sqrt(nearX * nearX + nearY * nearY);
      // The disc reads its own pixels: how far out a node is, and how far round.
      var stepR = seat.scale;
      var stepA = Math.max(1, budget) / (Math.PI * 2);
      var discPixels = {
        column: function (node) { return (radiusOf(view.x[node]) * stepR) | 0; },
        row: function (node) { return (angleOf(view.y[node]) * stepA) | 0; },
      };
      var picked = select(runs, budget, discPixels, view.shown, limits);
      view.shown = { child: picked.child, up: picked.up };
      view.picked = picked.count;

      ctx.strokeStyle = theme.branch;
      ctx.lineWidth = 1;
      ctx.beginPath();
      var drawn = 0;
      for (var each = 0; each < picked.count; each++) {
        var node = picked.child[each];
        var over = picked.up[each];
        var turn = angleOf(view.y[node]);
        var back = angleOf(view.y[over]);
        var out = radiusOf(view.x[node]) * seat.scale;
        var inward = radiusOf(view.x[over]) * seat.scale;
        var alongX = midX + Math.cos(turn) * inward;
        var alongY = midY + Math.sin(turn) * inward;
        // An arc, unless it is too flat to tell from a line. What decides that
        // is how far the arc bows away from its own chord, which for a small
        // turn is the radius times the square of it over eight: near the rim a
        // branch between neighbouring tips bows by a thousandth of a pixel, and
        // asking the canvas for an arc there costs nine times what a line does
        // and draws the same thing.
        if ((inward * (turn - back) * (turn - back)) / 8 < 0.5) {
          ctx.moveTo(midX + Math.cos(back) * inward, midY + Math.sin(back) * inward);
          ctx.lineTo(alongX, alongY);
        } else {
          ctx.moveTo(midX + Math.cos(back) * inward, midY + Math.sin(back) * inward);
          ctx.arc(midX, midY, inward, back, turn, turn < back);
        }
        ctx.moveTo(alongX, alongY);
        ctx.lineTo(midX + Math.cos(turn) * out, midY + Math.sin(turn) * out);
        drawn += 1;
      }
      ctx.stroke();

      var labels = 0;
      // Names, once a name has a whole row of pixels of rim to itself.
      var perTip = (Math.PI * 2 * seat.scale) / terminals();
      if (perTip >= LABEL_ROOM && picked.stride === 1) {
        ctx.fillStyle = theme.muted;
        ctx.font = Math.min(13, Math.max(9, perTip * 0.72)) + "px " + theme.font;
        ctx.textBaseline = "middle";
        for (var run = 0; run < runs.length; run++) {
          for (var i = runs[run][0]; i < runs[run][1]; i++) {
            var leaf = view.byRow[i];
            if (!view.length[leaf]) continue;
            var text = nameOf(leaf);
            if (!text) continue;
            var where = angleOf(view.y[leaf]);
            var edge = radiusOf(view.x[leaf]) * seat.scale + 4;
            var atX = midX + Math.cos(where) * edge;
            var atY = midY + Math.sin(where) * edge;
            if (atX < -40 || atX > box.wide + 40 || atY < -20 || atY > box.tall + 20) continue;
            ctx.textAlign = Math.cos(where) < 0 ? "right" : "left";
            ctx.fillText(text, atX, atY);
            labels += 1;
          }
        }
        ctx.textAlign = "left";
      }

      var cells = paintRings(ctx, theme, box, seat, midX, midY, runs);
      markFound(ctx, theme, box);
      var marked = box.wide >= RAIL_LEAST ? paintDial(ctx, theme, box, arc) : null;
      return {
        drawn: drawn,
        cells: cells,
        skipped: 0,
        labels: labels,
        stride: picked.stride,
        rowsInView: runs.reduce(function (sum, run) { return sum + run[1] - run[0]; }, 0),
        rail: marked,
      };
    }

    // The rootless walk, drawn. No arcs and no elbows: a branch is a straight
    // line between two points, because that is what the projection is.
    function paintSpread(theme) {
      var box = size();
      var ctx = canvas.getContext("2d");
      ctx.setTransform(box.ratio, 0, 0, box.ratio, 0, 0);
      ctx.clearRect(0, 0, box.wide, box.tall);

      var seat = discBox(box);
      var midX = seat.midX - disc.cx * seat.scale;
      var midY = seat.midY - disc.cy * seat.scale;
      var atX = function (node) { return midX + unitX(node) * seat.scale; };
      var atY = function (node) { return midY + unitY(node) * seat.scale; };

      // Nothing here maps a row to a pixel, so what a screen can tell apart is
      // its own edge: a tree drawn in a plane has its tips round the outside of
      // it, and there are as many places to put one as the window has pixels
      // round its rim.
      var budget = Math.max(256, Math.min(4000, Math.round((box.wide + box.tall) * 2)));
      var stepX = seat.scale;
      var pixel = {
        column: function (node) { return (unitX(node) * stepX) | 0; },
        row: function (node) { return (unitY(node) * stepX) | 0; },
      };
      var limits = {
        stop: null,
        // An edge with both ends off the same side of the window cannot cross
        // it. One with its ends off different sides might, and is kept.
        skip: function (walk, over) {
          var a = outside(atX(walk), atY(walk), box);
          var b = outside(atX(over), atY(over), box);
          return (a & b) !== 0;
        },
      };
      var picked = select([[0, view.byRow.length]], budget, pixel, view.shown, limits);
      view.shown = { child: picked.child, up: picked.up };
      view.picked = picked.count;

      ctx.strokeStyle = theme.branch;
      ctx.lineWidth = 1;
      ctx.beginPath();
      var drawn = 0;
      for (var each = 0; each < picked.count; each++) {
        var node = picked.child[each];
        var over = picked.up[each];
        ctx.moveTo(atX(over), atY(over));
        ctx.lineTo(atX(node), atY(node));
        drawn += 1;
      }
      ctx.stroke();

      var labels = 0;
      // Names, once the tips are far enough apart on screen to carry one.
      if (picked.stride === 1 && view.byRow.length * LABEL_ROOM < (box.wide + box.tall) * 2) {
        ctx.fillStyle = theme.muted;
        ctx.font = "11px " + theme.font;
        ctx.textBaseline = "middle";
        for (var i = 0; i < view.byRow.length; i++) {
          var leaf = view.byRow[i];
          if (!view.length[leaf]) continue;
          var text = nameOf(leaf);
          if (!text) continue;
          var px = atX(leaf), py = atY(leaf);
          if (outside(px, py, box)) continue;
          // Away from the middle, so a name leans out of the drawing rather
          // than across it.
          var lean = unitX(leaf) - (view.centreX || 0);
          ctx.textAlign = lean < 0 ? "right" : "left";
          ctx.fillText(text, px + (lean < 0 ? -4 : 4), py);
          labels += 1;
        }
        ctx.textAlign = "left";
      }

      markFound(ctx, theme, box);
      var marked = box.wide >= RAIL_LEAST ? paintDial(ctx, theme, box, null) : null;
      return {
        drawn: drawn,
        skipped: 0,
        labels: labels,
        stride: picked.stride,
        rowsInView: view.byRow.length,
        rail: marked,
      };
    }

    // ------------------------------------------------------------- the hand

    // How many terminals lie beyond each node, and the first and last row they
    // sit on. Worked out the first time something asks rather than at load,
    // because a tree nobody points at should not pay for it: on two million
    // nodes it is three arrays and one pass.
    //
    // The pass is bottom up by counting children off, not by walking, so it
    // does not care what order the nodes arrived in and cannot run out of
    // stack on a ladder.
    function beyond() {
      if (view.tipsBeyond) return;
      var count = view.count;
      var parent = view.parent;
      var kids = new Uint32Array(count);
      var at;
      for (at = 0; at < count; at++) {
        if (parent[at] !== 0xffffffff) kids[parent[at]] += 1;
      }
      var tips = new Uint32Array(count);
      var first = new Float32Array(count);
      var last = new Float32Array(count);
      var left = new Uint32Array(count);
      var queue = new Uint32Array(count);
      var head = 0, tail = 0;
      for (at = 0; at < count; at++) {
        left[at] = kids[at];
        first[at] = view.y[at];
        last[at] = view.y[at];
        if (!kids[at]) {
          tips[at] = 1;
          queue[tail++] = at;
        }
      }
      while (head < tail) {
        var node = queue[head++];
        var over = parent[node];
        if (over === 0xffffffff) continue;
        tips[over] += tips[node];
        if (first[node] < first[over]) first[over] = first[node];
        if (last[node] > last[over]) last[over] = last[node];
        left[over] -= 1;
        if (!left[over]) queue[tail++] = over;
      }
      view.tipsBeyond = tips;
      view.firstBeyond = first;
      view.lastBeyond = last;
    }

    // What is under a point, if anything is near enough to have been meant.
    //
    // The search is over what the last paint actually drew, which is both the
    // right answer and a small one: a reader can only point at what is on the
    // screen, and what is on the screen is a few thousand branches however many
    // the tree has.
    function at(px, py) {
      if (!view || !view.picked) return null;
      var box = size();
      var child = view.shown.child;
      var up = view.shown.up;
      var held = view.picked;
      var best = -1, near = 14 * 14;
      for (var scan = 0; scan < held; scan++) {
        for (var end = 0; end < 2; end++) {
          var node = end ? up[scan] : child[scan];
          var spot = whereOn(box, node);
          var dx = spot.x - px, dy = spot.y - py;
          var away = dx * dx + dy * dy;
          if (away < near) { near = away; best = node; }
        }
      }
      if (best < 0) return null;
      beyond();
      var spotted = whereOn(box, best);
      return {
        node: best,
        name: view.length[best] ? nameOf(best) : null,
        depth: view.x[best],
        tips: view.tipsBeyond[best],
        x: spotted.x,
        y: spotted.y,
      };
    }

    // Take the clade a node stands for: the rows its terminals sit on, or on
    // the disc the wedge they sweep. Without rows there is nothing to take, so
    // the rootless view goes to it instead.
    function focusOn(node) {
      if (!view || node === undefined || node === null) return false;
      beyond();
      var low = view.firstBeyond[node];
      var high = view.lastBeyond[node];
      var span = Math.max(3, (high - low) * 1.15);
      var middle = (low + high) / 2;
      if (mode === "rows") {
        camera.y0 = middle - span / 2;
        camera.y1 = middle + span / 2;
        return true;
      }
      var box = size();
      if (mode === "disc") {
        // The wedge, plus the room its own depth needs, so a clade fills the
        // window rather than sitting in a corner of it.
        var sweep = (span / Math.max(1, terminals())) * Math.PI * 2;
        var half = Math.max((Math.PI * 3) / terminals(), Math.min(1.1 * 1.5, sweep * 0.7));
        var turn = angleOf(middle);
        var out = (radiusOf(view.x[node]) + 1) / 2;
        disc.half = half;
        disc.cx = Math.cos(turn) * out;
        disc.cy = Math.sin(turn) * out;
        settleDisc();
        return true;
      }
      var spot = whereOn(box, node);
      var seat = discBox(box);
      disc.cx += (spot.x - seat.midX) / seat.scale;
      disc.cy += (spot.y - seat.midY) / seat.scale;
      disc.half = Math.max(disc.half / 3, (Math.PI * 3) / terminals());
      settleDisc();
      return true;
    }

    // ----------------------------------------------------------- the strips

    // A column of traits is drawn beside the names in the rectangular view and
    // as a ring outside the tips on the disc, which is what the crate does in
    // an SVG. The colours are not chosen here: they arrive with the layout,
    // resolved by the crate, so the strips on the canvas and the strips in the
    // exported figure cannot disagree about which blue is which.
    var STRIP_WIDE = 11;
    var STRIP_GAP = 2;

    function stripRoom() {
      var many = view && view.strips ? view.strips.length : 0;
      return many ? many * (STRIP_WIDE + STRIP_GAP) + STRIP_GAP : 0;
    }

    function stripInk(strip, node, dark) {
      var at = strip.of[node];
      if (at === 0xffffffff || at >= strip.levels.length) return null;
      var level = strip.levels[at];
      return dark ? level.dark : level.light;
    }

    // Beside the rows: one cell per column, at the row of every tip on screen.
    // A row too short to see is still drawn, because a strip is a count as much
    // as a picture and a gap in it says something.
    function paintStrips(ctx, theme, box, from, to, atY, perRow, wide) {
      var strips = view.strips;
      if (!strips || !strips.length) return 0;
      var tall = Math.max(1, perRow);
      var drawn = 0;
      for (var column = 0; column < strips.length; column++) {
        var strip = strips[column];
        var left = wide + STRIP_GAP + column * (STRIP_WIDE + STRIP_GAP);
        for (var at = from; at < to; at++) {
          var node = view.byRow[at];
          if (!view.length[node]) continue;
          var ink = stripInk(strip, node, theme.dark);
          if (!ink) continue;
          ctx.fillStyle = ink;
          ctx.fillRect(left, atY(view.y[node]) - tall / 2, STRIP_WIDE, tall);
          drawn += 1;
        }
      }
      return drawn;
    }

    // Round the rim: one ring per column, a wedge per tip. The same cells bent,
    // which is what the crate's trait rings are.
    function paintRings(ctx, theme, box, seat, midX, midY, runs) {
      var strips = view.strips;
      if (!strips || !strips.length) return 0;
      var many = terminals();
      var step = (Math.PI * 2) / Math.max(1, many);
      var drawn = 0;
      for (var column = 0; column < strips.length; column++) {
        var strip = strips[column];
        var inner = seat.scale * (1.02 + column * 0.05);
        var outer = seat.scale * (1.02 + column * 0.05 + 0.045);
        if (outer - inner < 1) continue;
        for (var run = 0; run < runs.length; run++) {
          for (var at = runs[run][0]; at < runs[run][1]; at++) {
            var node = view.byRow[at];
            if (!view.length[node]) continue;
            var ink = stripInk(strip, node, theme.dark);
            if (!ink) continue;
            var turn = angleOf(view.y[node]);
            ctx.fillStyle = ink;
            ctx.beginPath();
            ctx.arc(midX, midY, outer, turn - step / 2, turn + step / 2);
            ctx.arc(midX, midY, inner, turn + step / 2, turn - step / 2, true);
            ctx.closePath();
            ctx.fill();
            drawn += 1;
          }
        }
      }
      return drawn;
    }

    // ------------------------------------------------------------ selection

    // The branches to draw for a run of rows in a box `wide` by `tall`, as a
    // run of child and parent pairs.
    //
    // One row per pixel is all a screen can show, so past that they are stepped
    // over. Stepping alone draws a branch with nothing to hang it on: at any
    // zoom where the tree is taller than the canvas none of the kept rows is
    // the parent of another, and what is on screen is a hedge of loose
    // horizontal strokes rather than a tree. So every kept row is walked up to
    // the root and its ancestors are drawn with it.
    //
    // On most trees the walks meet almost at once and that is a handful of
    // extra branches. On a ladder it is not: every tip hangs off the spine, so
    // one walk is half the tree, and drawing a hundred and twenty thousand tip
    // caterpillar cost a hundred and fifty eight milliseconds a frame and did
    // not get cheaper when it was zoomed into. So the walk steps over an
    // ancestor that would be drawn inside the same pixel as the branch below
    // it, which is the same rule the rest of this crate uses along the other
    // axis: past one point per pixel the extra ones land on each other. The
    // ink is the same and the work is bounded by what a screen can hold.
    //
    // Everything drawn anywhere comes through here, which is what stops one
    // part of the canvas disagreeing with another about the shape of the tree.
    // `runs` is one or more [first, last) spans of rows. A circle can put the
    // rows in view either side of where its ends meet, and that is two spans of
    // one tree rather than two trees.
    function select(runs, budget, pixel, store, limits) {
      var whole = 0;
      for (var span = 0; span < runs.length; span++) whole += runs[span][1] - runs[span][0];
      var stride = Math.max(1, Math.ceil(whole / Math.max(1, budget)));
      view.visit += 1;
      var visit = view.visit;
      var seen = view.seen;
      var parent = view.parent;
      var child = store.child;
      var up = store.up;
      var count = 0;
      var sampled = 0;
      var column = pixel.column;
      var row = pixel.row;

      for (var pass = 0; pass < runs.length; pass++) {
      for (var at = runs[pass][0]; at < runs[pass][1]; at += stride) {
        var walk = view.byRow[at];
        sampled += 1;
        while (walk !== 0xffffffff && seen[walk] !== visit) {
          seen[walk] = visit;
          var over = parent[walk];
          if (over === 0xffffffff) break;
          while (
            parent[over] !== 0xffffffff &&
            seen[over] !== visit &&
            column(over) === column(walk) &&
            row(over) === row(walk)
          ) {
            seen[over] = visit;
            over = parent[over];
          }
          if (count === child.length) {
            var wideChild = new Uint32Array(child.length * 2);
            wideChild.set(child);
            child = wideChild;
            var wideUp = new Uint32Array(up.length * 2);
            wideUp.set(up);
            up = wideUp;
          }
          // A branch further out than anything the window can see is walked
          // past rather than drawn. The one that crosses into view is kept, so
          // the picture runs off the edge rather than stopping short of it.
          if (!limits || !limits.skip || !limits.skip(walk, over)) {
            child[count] = walk;
            up[count] = over;
            count += 1;
          }
          // And once a walk is further in than anything the window can see,
          // there is nothing left for it to draw. Every walk heads inward, so
          // this is where it ends.
          if (limits && limits.stop && limits.stop(over)) break;
          walk = over;
        }
      }
      }
      return { child: child, up: up, count: count, stride: stride, sampled: sampled };
    }

    // How the rectangular view sees its own pixels, worked out before the walk
    // rather than after it. The walk always reaches the root, so the left edge
    // is the root; and a parent is never deeper than its child, so the right
    // edge is the deepest of the rows sampled. Neither needs the walk to have
    // run.
    function rowPixels(runs, wide, tall, spanY, stride) {
      var x = view.x;
      var y = view.y;
      var lowX = x[view.root];
      var highX = lowX;
      for (var run = 0; run < runs.length; run++) {
        for (var look = runs[run][0]; look < runs[run][1]; look += stride) {
          var seenX = x[view.byRow[look]];
          if (seenX > highX) highX = seenX;
        }
      }
      var acrossX = Math.max(1, wide) / (highX - lowX || 1);
      var downY = Math.max(1, tall) / (spanY || 1);
      return {
        column: function (node) { return ((x[node] - lowX) * acrossX) | 0; },
        row: function (node) { return (y[node] * downY) | 0; },
      };
    }

    // The span in x that holds a selection, which is the root on the left and
    // the deepest branch in it on the right. Both endpoints of every branch,
    // because a parent stepped over by the pixel test is still drawn to.
    function spread(picked) {
      var lowX = Infinity, highX = -Infinity;
      for (var scan = 0; scan < picked.count; scan++) {
        var a = view.x[picked.child[scan]];
        var c = view.x[picked.up[scan]];
        if (a < lowX) lowX = a;
        if (c < lowX) lowX = c;
        if (a > highX) highX = a;
        if (c > highX) highX = c;
      }
      if (!(highX > lowX)) return { lowX: view.bounds.lowX, highX: view.bounds.highX };
      return { lowX: lowX, highX: highX };
    }

    function paint(theme) {
      if (!view) return { drawn: 0, skipped: 0 };
      if (mode === "disc") return paintDisc(theme);
      if (mode === "spread") return paintSpread(theme);
      var box = size();
      var ctx = canvas.getContext("2d");
      ctx.setTransform(box.ratio, 0, 0, box.ratio, 0, 0);
      ctx.clearRect(0, 0, box.wide, box.tall);

      var strip = rail(box);
      // The tree draws into what is left when the rail and the trait columns
      // have taken theirs.
      var wide = (strip ? strip.x0 : box.wide) - stripRoom();

      var spanX = camera.x1 - camera.x0 || 1;
      var spanY = camera.y1 - camera.y0 || 1;
      var sx = wide / spanX;
      var sy = box.tall / spanY;
      var atX = function (value) { return (value - camera.x0) * sx; };
      var atY = function (value) { return (value - camera.y0) * sy; };

      var from = firstRow(camera.y0 - 1);
      var to = firstRow(camera.y1 + 1);
      var wanted = to - from;
      var runs = [[from, to]];
      var stride = Math.max(1, Math.ceil((to - from) / Math.max(1, box.tall)));
      var picked = select(
        runs,
        box.tall,
        rowPixels(runs, wide, box.tall, spanY, stride),
        view.shown
      );
      view.shown = { child: picked.child, up: picked.up };
      var count = picked.count;
      stride = picked.stride;

      // The depth axis follows what is drawn rather than being zoomed alongside
      // the rows. Zooming both narrowed the window in x as well until it held
      // no branches at all: past fifty wheel notches on a million tip tree the
      // canvas went blank while the rows were still there. Since the walk above
      // always reaches the root, the left edge is the root and the right edge
      // is the deepest tip on screen, so the axis stands still while it is
      // panned and only gives ground back as a clade is entered.
      var reach = spread(picked);
      var margin = (reach.highX - reach.lowX) * 0.28;
      camera.x0 = reach.lowX - margin * 0.05;
      camera.x1 = reach.highX + margin;
      spanX = camera.x1 - camera.x0 || 1;
      sx = wide / spanX;

      view.picked = count;

      ctx.strokeStyle = theme.branch;
      ctx.lineWidth = 1;
      ctx.beginPath();
      var drawn = 0;
      for (var each = 0; each < count; each++) {
        var node = picked.child[each];
        var over = picked.up[each];
        var y = atY(view.y[node]);
        var x1 = atX(view.x[node]);
        var x0 = atX(view.x[over]);
        ctx.moveTo(x0, y);
        ctx.lineTo(x1, y);
        // The elbow up to the parent's own row.
        ctx.moveTo(x0, y);
        ctx.lineTo(x0, atY(view.y[over]));
        drawn += 1;
      }
      ctx.stroke();

      // Names, once there is room for them. Below that a name would be drawn on
      // top of its neighbour and read as neither.
      var perRow = sy;
      var labels = 0;
      if (perRow >= LABEL_ROOM && stride === 1) {
        ctx.fillStyle = theme.muted;
        ctx.font = Math.min(13, Math.max(9, perRow * 0.72)) + "px " + theme.font;
        ctx.textBaseline = "middle";
        // Names stop where the rail begins. Without the stop they run under it
        // and the silhouette is drawn over their tails, which reads as a name
        // that has been cut rather than as a name behind something.
        for (var i = from; i < to; i++) {
          var leaf = view.byRow[i];
          if (!view.length[leaf]) continue;
          var text = nameOf(leaf);
          if (!text) continue;
          var left = atX(view.x[leaf]) + 4;
          if (left + ctx.measureText(text).width > wide) continue;
          ctx.fillText(text, left, atY(view.y[leaf]));
          labels += 1;
        }
      }
      var cells = paintStrips(ctx, theme, box, from, to, atY, sy, wide);
      markFound(ctx, theme, box);
      var marked = strip ? paintRail(ctx, theme, box, strip) : null;

      return {
        drawn: drawn,
        cells: cells,
        skipped: wanted - picked.sampled,
        labels: labels,
        stride: stride,
        rowsInView: wanted,
        rail: marked,
      };
    }

    // The rail: the whole tree at the height of the canvas, with the rows on
    // screen marked on it. Its own picture never moves, so what a reader
    // follows is the mark travelling down a shape that stays put.
    function paintRail(ctx, theme, box, strip) {
      var seen = overview(box, theme);

      // The edge it stands behind, so the rail reads as a margin and not as
      // more tree.
      ctx.strokeStyle = theme.frame;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(strip.x0 + 0.5, 0);
      ctx.lineTo(strip.x0 + 0.5, box.tall);
      ctx.stroke();

      if (seen.plate) ctx.drawImage(seen.plate, strip.x0, 0, strip.wide, box.tall);
      else strokeOverview(ctx, theme, box, strip, seen);

      // The rows on screen. Three rows out of two million is a thousandth of a
      // pixel, so the mark is held to a size a reader can see and a hand can
      // catch, and it is kept on the rail rather than allowed to hang off it.
      var top = heightOfRow(camera.y0, box);
      var foot = heightOfRow(camera.y1, box);
      var deep = Math.max(MARK_LEAST, foot - top);
      if (top + deep > box.tall) top = box.tall - deep;
      if (top < 0) top = 0;
      ctx.fillStyle = theme.window;
      ctx.fillRect(strip.x0 + 1, top, strip.wide - 1, deep);
      ctx.strokeStyle = theme.edge;
      ctx.strokeRect(strip.x0 + 1.5, top + 0.5, strip.wide - 2, Math.max(1, deep - 1));

      return { x0: strip.x0, wide: strip.wide, top: top, deep: deep, drawn: seen.count };
    }

    // --------------------------------------------------------------- camera

    // Only the rows. The depth follows them, in `paint`.
    function zoomAt(px, py, factor) {
      var box = size();
      if (mode !== "rows") {
        // Keep whatever is under the hand under the hand, which is what a map
        // does and what a run of rows cannot do because it has only one axis.
        var seat = discBox(box);
        var atX = disc.cx + (px - seat.midX) / seat.scale;
        var atY = disc.cy + (py - seat.midY) / seat.scale;
        // A map keeps whatever is under the hand under the hand. The middle of
        // a disc is a hole with nothing in it, and it is exactly where a hand
        // rests when the whole circle is on screen, so a wheel there converged
        // on emptiness. Where there is nothing to hold on to, hold on to the
        // nearest thing there is, which is the ring the root sits on.
        if (mode === "disc") {
          var out = Math.sqrt(atX * atX + atY * atY);
          if (out < DISC_HOLE) {
            if (out > 1e-9) {
              atX = (atX / out) * DISC_HOLE;
              atY = (atY / out) * DISC_HOLE;
            } else {
              // Dead on the middle there is no direction to push out along, so
              // it goes to where the root sits. Any other bearing is a guess,
              // and on a tree whose root is round the other side it is a guess
              // at a piece of the ring with nothing on it.
              var facing = angleOf(view.y[view.root]);
              atX = Math.cos(facing) * DISC_HOLE;
              atY = Math.sin(facing) * DISC_HOLE;
            }
          }
        }
        // The same two stops the rows have, read in the other geometry. Out as
        // far as the whole drawing and half again, and in as far as three tips
        // of rim, which is the circular reading of the three rows a run of them
        // stops at. Without the second one the wheel went on going in for ever,
        // long past anything being left to see.
        var half = disc.half / factor;
        var widest = 1.1 * 1.5;
        var closest = Math.max(1e-6, (Math.PI * 3) / terminals());
        if (half > widest) half = widest;
        if (half < closest) half = closest;
        var wasHalf = disc.half, wasX = disc.cx, wasY = disc.cy;
        // Kept so a frame that comes back with nothing on it can be undone.
        // Whether a window holds any ink is not a question its outline can
        // answer: a window can sit squarely on the drawing and land in the gap
        // between two branches. The drawing itself knows, so the drawing is
        // what decides, one frame later.
        view.wasCamera = { cx: wasX, cy: wasY, half: wasHalf };
        disc.half = half;
        var after = discBox(box);
        disc.cx = atX - (px - after.midX) / after.scale;
        disc.cy = atY - (py - after.midY) / after.scale;
        settleDisc();
        // And a wheel that takes the window off the drawing altogether has gone
        // too far. The anchor above is what keeps the disc's middle from being
        // a hole to fall into; this catches the other way out, which is zooming
        // at the very edge until the drawing is behind you.
        if (!anythingInView(box)) {
          disc.half = wasHalf;
          disc.cx = wasX;
          disc.cy = wasY;
          return false;
        }
        return true;
      }
      var fy = py / box.tall;
      var atY = camera.y0 + (camera.y1 - camera.y0) * fy;
      var tall = (camera.y1 - camera.y0) / factor;
      var b = view.bounds;
      // Out as far as the whole tree and half again, in as far as three rows,
      // which is a cherry and its parent.
      var tallest = (b.highY - b.lowY + 1) * 1.5;
      if (tall > tallest) tall = tallest;
      if (tall < 3) tall = 3;
      camera.y0 = atY - tall * fy;
      camera.y1 = atY + tall * (1 - fy);
      return true;
    }

    // Half a screen of empty is as far as the rows go. Without a stop the tree
    // can be pushed off the canvas altogether, and white with no rows on it
    // gives a reader nothing to drag back by.
    function settle(y0, span) {
      var b = view.bounds;
      var first = b.lowY - span * 0.5;
      var last = b.highY - span * 0.5;
      if (y0 < first) y0 = first;
      if (y0 > last) y0 = last;
      camera.y0 = y0;
      camera.y1 = y0 + span;
    }

    // Rows only. A drag sideways is taken and dropped: the depth axis is fitted
    // to what is drawn on every paint, so nudging it would be undone before the
    // frame was on screen.
    function panBy(dx, dy) {
      var box = size();
      if (mode !== "rows") {
        var seat = discBox(box);
        disc.cx -= dx / seat.scale;
        disc.cy -= dy / seat.scale;
        settleDisc();
        return;
      }
      var span = camera.y1 - camera.y0;
      settle(camera.y0 - (dy / box.tall) * span, span);
    }

    // Whether the window still overlaps the drawing at all. This is a test of
    // the outline and not of the ink: on the disc the ink lies between the hole
    // and the rim, so a window wholly inside the one or wholly outside the
    // other holds nothing for certain, and the rootless walk fills a box a
    // window can miss entirely. A window that passes can still land in a gap
    // between two branches, which is what any map does and what the projection
    // itself puts there.
    function anythingInView(box) {
      var seat = discBox(box);
      var x0 = disc.cx - seat.halfW, x1 = disc.cx + seat.halfW;
      var y0 = disc.cy - disc.half, y1 = disc.cy + disc.half;
      if (mode === "disc") {
        // Cheapest first: the ink lies between the hole and the rim, so a
        // window wholly inside the one or wholly outside the other holds
        // nothing and there is no need to look further.
        var nearX = Math.max(x0, Math.min(0, x1));
        var nearY = Math.max(y0, Math.min(0, y1));
        var farX = Math.max(Math.abs(x0), Math.abs(x1));
        var farY = Math.max(Math.abs(y0), Math.abs(y1));
        var nearest = Math.sqrt(nearX * nearX + nearY * nearY);
        var farthest = Math.sqrt(farX * farX + farY * farY);
        if (!(farthest >= DISC_HOLE && nearest <= 1)) return false;
      } else if (!(x0 <= 1 && x1 >= -1 && y0 <= 1 && y1 >= -1)) {
        return false;
      }
      return true;
    }

    // The disc can be pushed until the rim is off the canvas, but not until
    // there is nothing left to aim at.
    function settleDisc() {
      if (disc.cx < -1.2) disc.cx = -1.2;
      if (disc.cx > 1.2) disc.cx = 1.2;
      if (disc.cy < -1.2) disc.cy = -1.2;
      if (disc.cy > 1.2) disc.cy = 1.2;
    }

    // Is this point, in canvas pixels, on the small picture of the whole tree
    // rather than on the tree itself? The rail in rows, the dial on the disc.
    function onMap(px, py) {
      if (!view) return false;
      var box = size();
      if (mode !== "rows") {
        var spot = dial(box);
        return (
          !!spot &&
          px >= spot.x0 && px <= spot.x0 + spot.side &&
          py >= spot.y0 && py <= spot.y0 + spot.side
        );
      }
      var strip = rail(box);
      return !!strip && px >= strip.x0;
    }

    // Put the rows on screen where this height on the rail points, keeping how
    // many of them there are. A click and a drag are the same gesture: the
    // window follows the hand rather than being nudged by it.
    function jumpTo(px, py) {
      if (!view) return;
      var box = size();
      if (mode !== "rows") {
        var spot = dial(box);
        var held = view.dial;
        if (!spot || !held) return;
        disc.cx = (px - (spot.x0 + spot.side / 2)) / held.span;
        disc.cy = (py - (spot.y0 + spot.side / 2)) / held.span;
        settleDisc();
        return;
      }
      var span = camera.y1 - camera.y0;
      settle(rowAtHeight(py, box) - span / 2, span);
    }

    // Puts a named tip in the middle, without changing how much is on screen.
    // Take the view to a tip by name, in whichever projection is being looked
    // through. It used to move the row camera and nothing else, so in the two
    // projections that are driven by the other camera a hit moved nothing at
    // all and still reported success.
    // Every tip whose name answers to `text`, as node indices in row order.
    //
    // Over the bytes rather than over strings: the names arrive as one blob and
    // decoding each of a million of them to compare it took nearly two seconds,
    // which is a search a reader gives up on. Comparing bytes allocates
    // nothing, and folding ASCII case as it goes costs one branch.
    //
    // `how` is "exact", "loose" for case insensitive, "starts" for a prefix, or
    // "in" for anywhere in the name. A list of names separated by commas,
    // spaces or newlines is taken as a list and every one of them looked for.
    function find(text, how) {
      if (!view) return new Uint32Array(0);
      var wanted = String(text || "").trim();
      if (!wanted) return new Uint32Array(0);
      var many = wanted.split(/[\s,]+/).filter(Boolean);
      if (many.length > 1) {
        var seen = new Uint8Array(view.count);
        for (var each = 0; each < many.length; each++) {
          var one = find(many[each], how);
          for (var mark = 0; mark < one.length; mark++) seen[one[mark]] = 1;
        }
        var all = [];
        for (var walk = 0; walk < view.byRow.length; walk++) {
          if (seen[view.byRow[walk]]) all.push(view.byRow[walk]);
        }
        return Uint32Array.from(all);
      }

      var fold = how !== "exact";
      var query = encoder.encode(fold ? wanted.toLowerCase() : wanted);
      var names = view.names;
      var start = view.start;
      var length = view.length;
      var hits = [];
      var lower = function (byte) {
        return byte >= 65 && byte <= 90 ? byte + 32 : byte;
      };
      for (var order = 0; order < view.byRow.length; order++) {
        var node = view.byRow[order];
        var len = length[node];
        if (!len) continue;
        var at = start[node];
        if (how === "exact" || how === "loose") {
          if (len !== query.length) continue;
          var same = true;
          for (var i = 0; i < len && same; i++) {
            var here = fold ? lower(names[at + i]) : names[at + i];
            if (here !== query[i]) same = false;
          }
          if (same) hits.push(node);
          continue;
        }
        if (len < query.length) continue;
        var last = how === "starts" ? 0 : len - query.length;
        for (var from = 0; from <= last; from++) {
          var run = true;
          for (var step = 0; step < query.length && run; step++) {
            if (lower(names[at + from + step]) !== query[step]) run = false;
          }
          if (run) { hits.push(node); break; }
        }
      }
      return Uint32Array.from(hits);
    }

    // Take the view to a tip by name, in whichever projection is being looked
    // through. It used to move the row camera and nothing else, so in the two
    // projections that are driven by the other camera a hit moved nothing at
    // all and still reported success.
    function goTo(name, how) {
      var hits = find(name, how || "exact");
      if (!hits.length) return false;
      view.found = hits;
      view.at = 0;
      settleOn(hits[0]);
      return true;
    }

    // Move to the next of what the last search found, coming round at the end.
    function nextFound(step) {
      if (!view || !view.found || !view.found.length) return null;
      var many = view.found.length;
      view.at = ((view.at + (step || 1)) % many + many) % many;
      settleOn(view.found[view.at]);
      return { at: view.at, of: many };
    }

    function settleOn(node) {
      if (mode === "rows") {
        // Close enough that the tip has a name beside it, rather than
        // wherever the reader happened to be zoomed to.
        var tall = Math.min(camera.y1 - camera.y0, 60);
        camera.y0 = view.y[node] - tall * 0.5;
        camera.y1 = view.y[node] + tall * 0.5;
        return;
      }
      disc.cx = unitX(node);
      disc.cy = unitY(node);
      if (disc.half > 0.08) disc.half = 0.08;
      settleDisc();
    }

    // Where the last search landed, so a paint can put a ring round it. A tip
    // found and not marked is a tip the reader still has to hunt for.
    function markFound(ctx, theme, box) {
      var hits = view.found;
      if (!hits || !hits.length) return;
      ctx.strokeStyle = theme.edge;
      ctx.lineWidth = 2;
      var drawn = 0;
      for (var each = 0; each < hits.length && drawn < 400; each++) {
        var at = whereOn(box, hits[each]);
        if (!at) continue;
        if (at.x < -20 || at.x > box.wide + 20 || at.y < -20 || at.y > box.tall + 20) continue;
        ctx.beginPath();
        ctx.arc(at.x, at.y, each === view.at ? 8 : 5, 0, Math.PI * 2);
        ctx.stroke();
        drawn += 1;
      }
      ctx.lineWidth = 1;
    }

    // Where a node lands, given a box already measured.
    function whereOn(box, node) {
      if (mode !== "rows") {
        var seat = discBox(box);
        return {
          x: seat.midX + (unitX(node) - disc.cx) * seat.scale,
          y: seat.midY + (unitY(node) - disc.cy) * seat.scale,
        };
      }
      // The same width the paint drew into, trait columns taken off it as
      // well: without that the hand would name a branch a little to the right
      // of the one under it.
      var strip = rail(box);
      var wide = (strip ? strip.x0 : box.wide) - stripRoom();
      return {
        x: ((view.x[node] - camera.x0) / (camera.x1 - camera.x0 || 1)) * wide,
        y: ((view.y[node] - camera.y0) / (camera.y1 - camera.y0 || 1)) * box.tall,
      };
    }

    // What the window is looking at, in rows, which is what a figure of this
    // view would be asked for.
    function looking() {
      if (!view) return null;
      var from, to;
      if (mode !== "rows") {
        var runs = mode === "disc" ? runsForArc(arcInView(size())) : [[0, view.byRow.length]];
        // Two runs mean the window straddles where the circle's ends meet, and
        // the pair a reader wants named is the whole of what is on screen.
        from = runs[0][0];
        to = runs[runs.length - 1][1];
        if (runs.length > 1) { from = 0; to = view.byRow.length; }
      } else {
        from = firstRow(camera.y0);
        to = firstRow(camera.y1);
      }
      // In from each end rather than along the whole run. Walking every row in
      // view and decoding its name is fine on a screenful and is not fine on
      // the first sight of a million tip tree, where the whole file is in view:
      // it took 3.3 seconds, and it was doing it on every wheel notch.
      var first = null, last = null;
      for (var at = from; at < to && first === null; at++) {
        var node = view.byRow[at];
        if (view.length[node]) first = nameOf(node);
      }
      for (var back = to - 1; back >= from && last === null; back--) {
        var end = view.byRow[back];
        if (view.length[end]) last = nameOf(end);
      }
      // How many rows the window holds, not how many nodes fall inside them.
      // A parent sits on a row of its own between its children, so the node
      // count is about twice the row count and saying it was rows was wrong.
      // How many rows the tree has. With a root that is the span of the row
      // numbers; without one there are no row numbers at all, because `bounds`
      // there holds coordinates, and the count is simply the terminals. Reading
      // the coordinates as rows said a tree of twenty thousand tips had one,
      // and the figure it asked for was folded to eight.
      var b = view.bounds;
      var whole = view.rootless ? view.byRow.length : b.highY - b.lowY + 1;
      var span;
      if (mode === "rows") {
        span = Math.min(camera.y1, b.highY + 0.5) - Math.max(camera.y0, b.lowY - 0.5);
      } else if (mode === "disc") {
        // A wedge of the circle is a run of rows, so a reader zoomed into one
        // is not looking at the whole tree and the figure they ask for should
        // not be drawn as though they were.
        var seen = arcInView(size());
        span = seen ? (seen.span / (Math.PI * 2)) * whole : whole;
      } else {
        // The rootless walk holds every tip at every zoom: there is no run of
        // rows to have only some of.
        span = whole;
      }
      return { first: first, last: last, rows: Math.max(0, Math.round(span)) };
    }

    return {
      load: load,
      paint: paint,
      home: home,
      zoomAt: zoomAt,
      panBy: panBy,
      onMap: onMap,
      jumpTo: jumpTo,
      // Undo the last wheel notch, for a caller that has just painted with it
      // and found nothing there. It can only be used once per notch, so a
      // caller cannot walk itself backwards for ever.
      stepBack: function () {
        if (!view || mode === "rows" || !view.wasCamera) return false;
        disc.cx = view.wasCamera.cx;
        disc.cy = view.wasCamera.cy;
        disc.half = view.wasCamera.half;
        view.wasCamera = null;
        return true;
      },
      // Which projection is being looked through: "rows", "disc" or "spread".
      // The first two are the same layout read two ways and cost nothing to
      // change between. The third is a different walk and the page has to have
      // asked the program for it, so it is refused until the coordinates for
      // it have arrived.
      shape: function (want) {
        if (!view || want === mode) return mode;
        if (want === "spread" && !view.rootless) return mode;
        if (want !== "spread" && view.rootless) return mode;
        if (want !== "rows" && want !== "disc" && want !== "spread") return mode;
        mode = want;
        home();
        return mode;
      },
      shapeNow: function () { return mode; },
      // The columns as they were resolved, for a key beside the picture.
      strips: function () {
        return view && view.strips ? view.strips : [];
      },
      rootless: function () { return !!(view && view.rootless); },
      // Where a node sits on the canvas, in the projection being looked
      // through. The one place that answers it, so a check of a projection can
      // ask whether the canvas agrees with the crate's own arithmetic.
      where: function (node) {
        if (!view) return null;
        return whereOn(size(), node);
      },
      goTo: goTo,
      find: find,
      nextFound: nextFound,
      found: function () {
        if (!view || !view.found || !view.found.length) return { count: 0, at: 0, names: [] };
        // The names, for a figure to be drawn with them marked. Capped, because
        // a command line is a thing a person reads and a search can match a
        // hundred thousand tips.
        var names = [];
        for (var each = 0; each < view.found.length && names.length < 12; each++) {
          if (view.length[view.found[each]]) names.push(nameOf(view.found[each]));
        }
        return { count: view.found.length, at: view.at || 0, names: names };
      },
      at: at,
      focusOn: focusOn,
      looking: looking,
      loaded: function () { return !!view; },
      count: function () { return view ? view.count : 0; },
      // What the last paint put on the canvas, and the depth it fitted them
      // into. Both are answers about the picture rather than about the tree,
      // which is what a check of the picture needs to ask.
      // The branches the last paint put on the canvas, as child and parent
      // pairs. An answer about the picture rather than about the tree, which is
      // what a check of the picture needs to ask.
      shown: function () {
        if (!view) return { child: new Uint32Array(0), up: new Uint32Array(0), count: 0 };
        var held = view.picked || 0;
        return {
          child: view.shown.child.subarray(0, held),
          up: view.shown.up.subarray(0, held),
          count: held,
        };
      },
      // How far across the drawing the window reaches, in whatever the
      // projection measures across in. There is no row camera in the two flat
      // views, so asking for one there used to throw.
      depth: function () {
        if (!view) return null;
        if (mode !== "rows") {
          var seat = discBox(size());
          return { x0: disc.cx - seat.halfW, x1: disc.cx + seat.halfW };
        }
        return { x0: camera.x0, x1: camera.x1 };
      },
    };
  }

  return { make: make };
})();
