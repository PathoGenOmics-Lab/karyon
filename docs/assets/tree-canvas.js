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
  var RAIL_PAD = 5;
  var RAIL_LEAST = 420;
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

    function load(placed) {
      view = placed;
      // A layout with no root sends the order to read its terminals in, since
      // it has no rows for the page to sort by. That order is what stands in
      // for rows here: it is what one is picked out of every few from.
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
      return { x0: box.wide - RAIL_WIDE, wide: RAIL_WIDE, tall: box.tall };
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
      if (view.overviewFor === box.tall && view.overviewInk === theme.faint) {
        return view.overview;
      }
      var runs = [[0, view.byRow.length]];
      var stride = Math.max(1, Math.ceil(view.byRow.length / Math.max(1, box.tall)));
      var picked = select(
        runs,
        box.tall,
        rowPixels(runs, RAIL_WIDE, box.tall, view.bounds.highY - view.bounds.lowY + 1, stride),
        { child: new Uint32Array(1 << 13), up: new Uint32Array(1 << 13) }
      );
      var seen = {
        child: picked.child,
        up: picked.up,
        count: picked.count,
        reach: spread(picked),
      };
      seen.plate = plate(box, theme, seen);
      view.overview = seen;
      view.overviewFor = box.tall;
      view.overviewInk = theme.faint;
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
      sheet.width = Math.max(1, Math.round(RAIL_WIDE * ratio));
      sheet.height = Math.max(1, Math.round(box.tall * ratio));
      var ink = sheet.getContext("2d");
      if (!ink) return null;
      ink.setTransform(ratio, 0, 0, ratio, 0, 0);
      strokeOverview(ink, theme, box, { x0: 0, wide: RAIL_WIDE }, seen);
      return sheet;
    }

    // The silhouette itself, given somewhere to put it.
    function strokeOverview(ink, theme, box, strip, seen) {
      var inner = RAIL_WIDE - RAIL_PAD * 2;
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
    function dialPlate(theme) {
      var papers = canvas.ownerDocument;
      var span = DIAL_SIDE / 2 / 1.1;
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
      sheet.width = Math.round(DIAL_SIDE * ratio);
      sheet.height = Math.round(DIAL_SIDE * ratio);
      var ink = sheet.getContext("2d");
      if (!ink) return held;
      ink.setTransform(ratio, 0, 0, ratio, 0, 0);
      strokeDisc(ink, theme, held, DIAL_SIDE / 2, DIAL_SIDE / 2, span);
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
      return {
        x0: box.wide - DIAL_SIDE - DIAL_EDGE,
        y0: box.tall - DIAL_SIDE - DIAL_EDGE,
        side: DIAL_SIDE,
      };
    }

    function paintDial(ctx, theme, box, arc) {
      var spot = dial(box);
      if (!spot) return null;
      if (!view.dial || view.dialInk !== theme.faint || view.dialMode !== mode) {
        view.dial = dialPlate(theme);
        view.dialInk = theme.faint;
        view.dialMode = mode;
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

      var marked = box.wide >= RAIL_LEAST ? paintDial(ctx, theme, box, arc) : null;
      return {
        drawn: drawn,
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
      // The tree draws into what is left when the rail has taken its width.
      var wide = strip ? strip.x0 : box.wide;

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
      var marked = strip ? paintRail(ctx, theme, box, strip) : null;

      return {
        drawn: drawn,
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
        var half = disc.half / factor;
        if (half > 1.4) half = 1.4;
        if (half < 1e-5) half = 1e-5;
        disc.half = half;
        var after = discBox(box);
        disc.cx = atX - (px - after.midX) / after.scale;
        disc.cy = atY - (py - after.midY) / after.scale;
        settleDisc();
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
    function goTo(name) {
      if (!view) return false;
      for (var node = 0; node < view.count; node++) {
        if (!view.length[node]) continue;
        if (nameOf(node) !== name) continue;
        var tall = camera.y1 - camera.y0;
        camera.y0 = view.y[node] - tall * 0.5;
        camera.y1 = view.y[node] + tall * 0.5;
        return true;
      }
      return false;
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
      return { first: first, last: last, rows: to - from };
    }

    return {
      load: load,
      paint: paint,
      home: home,
      zoomAt: zoomAt,
      panBy: panBy,
      onMap: onMap,
      jumpTo: jumpTo,
      // Which projection is being looked through. The layout does not change:
      // the same rows and depths become angles and radii.
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
      rootless: function () { return !!(view && view.rootless); },
      // Where a node sits on the canvas, in the projection being looked
      // through. The one place that answers it, so a check of the disc can ask
      // whether the canvas agrees with the crate's own arithmetic.
      where: function (node) {
        if (!view) return null;
        var box = size();
        if (mode !== "rows") {
          // Both flat projections put a node somewhere in the same two units,
          // so there is one answer here and not two.
          var seat = discBox(box);
          return {
            x: seat.midX + (unitX(node) - disc.cx) * seat.scale,
            y: seat.midY + (unitY(node) - disc.cy) * seat.scale,
          };
        }
        var strip = rail(box);
        var wide = strip ? strip.x0 : box.wide;
        return {
          x: ((view.x[node] - camera.x0) / (camera.x1 - camera.x0 || 1)) * wide,
          y: ((view.y[node] - camera.y0) / (camera.y1 - camera.y0 || 1)) * box.tall,
        };
      },
      goTo: goTo,
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
      depth: function () { return { x0: camera.x0, x1: camera.x1 }; },
    };
  }

  return { make: make };
})();
