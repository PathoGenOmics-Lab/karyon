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

  function make(canvas) {
    var view = null;
    var camera = null;
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
      view.byRow = order(placed);
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
      view.seen = new Int32Array(placed.count);
      view.visit = 0;
      view.shown = new Uint32Array(1 << 14);
      home();
    }

    function home() {
      var b = view.bounds;
      // A margin on the right for the names, which are drawn outward from the
      // tip and are not in the coordinates.
      camera = { x0: b.lowX, x1: b.highX + (b.highX - b.lowX) * 0.25, y0: b.lowY - 0.5, y1: b.highY + 0.5 };
    }

    function size() {
      var ratio = window.devicePixelRatio || 1;
      var wide = canvas.clientWidth || 1;
      var tall = canvas.clientHeight || 1;
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
    // only if the canvas changes height, because it is the same picture at a
    // different resolution and not a second opinion about the tree: it comes
    // out of the same `select` the main view is drawn from.
    function overview(box) {
      if (view.overviewFor === box.tall) return view.overview;
      var picked = select(0, view.byRow.length, box.tall, new Uint32Array(1 << 13));
      view.overview = { nodes: picked.nodes, count: picked.count, reach: spread(picked.nodes, picked.count) };
      view.overviewFor = box.tall;
      return view.overview;
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

    // ------------------------------------------------------------ selection

    // The nodes to draw for a run of rows in `tall` pixels of height, as a run
    // of indices into `into`.
    //
    // One row per pixel is all a screen can show, so past that they are stepped
    // over. Stepping alone draws a branch with nothing to hang it on: at any
    // zoom where the tree is taller than the canvas none of the kept rows is
    // the parent of another, and what is on screen is a hedge of loose
    // horizontal strokes rather than a tree. So every kept row is walked up to
    // the root and its ancestors are drawn with it. The walks meet and stop, so
    // it is a handful of extra nodes, and it is what puts a trunk on the
    // picture.
    //
    // Everything drawn anywhere comes through here, which is what stops one
    // part of the canvas disagreeing with another about the shape of the tree.
    function select(from, to, tall, into) {
      var stride = Math.max(1, Math.ceil((to - from) / Math.max(1, tall)));
      view.visit += 1;
      var visit = view.visit;
      var seen = view.seen;
      var shown = into;
      var count = 0;
      var sampled = 0;
      for (var at = from; at < to; at += stride) {
        var walk = view.byRow[at];
        sampled += 1;
        while (walk !== 0xffffffff && seen[walk] !== visit) {
          seen[walk] = visit;
          if (count === shown.length) {
            var wider = new Uint32Array(shown.length * 2);
            wider.set(shown);
            shown = wider;
          }
          shown[count] = walk;
          count += 1;
          walk = view.parent[walk];
        }
      }
      return { nodes: shown, count: count, stride: stride, sampled: sampled };
    }

    // The span in x that holds a selection, which is the root on the left and
    // the deepest node in it on the right.
    function spread(nodes, count) {
      var lowX = Infinity, highX = -Infinity;
      for (var scan = 0; scan < count; scan++) {
        var it = nodes[scan];
        if (view.x[it] < lowX) lowX = view.x[it];
        if (view.x[it] > highX) highX = view.x[it];
      }
      if (!(highX > lowX)) return { lowX: view.bounds.lowX, highX: view.bounds.highX };
      return { lowX: lowX, highX: highX };
    }

    function paint(theme) {
      if (!view) return { drawn: 0, skipped: 0 };
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
      var picked = select(from, to, box.tall, view.shown);
      view.shown = picked.nodes;
      var shown = picked.nodes;
      var count = picked.count;
      var stride = picked.stride;

      // The depth axis follows what is drawn rather than being zoomed alongside
      // the rows. Zooming both narrowed the window in x as well until it held
      // no branches at all: past fifty wheel notches on a million tip tree the
      // canvas went blank while the rows were still there. Since the walk above
      // always reaches the root, the left edge is the root and the right edge
      // is the deepest tip on screen, so the axis stands still while it is
      // panned and only gives ground back as a clade is entered.
      var reach = spread(shown, count);
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
        var node = shown[each];
        var up = view.parent[node];
        if (up === 0xffffffff) continue;
        var y = atY(view.y[node]);
        var x1 = atX(view.x[node]);
        var x0 = atX(view.x[up]);
        ctx.moveTo(x0, y);
        ctx.lineTo(x1, y);
        // The elbow up to the parent's own row.
        var py = atY(view.y[up]);
        ctx.moveTo(x0, y);
        ctx.lineTo(x0, py);
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
        for (var i = from; i < to; i++) {
          var leaf = view.byRow[i];
          if (!view.length[leaf]) continue;
          var text = nameOf(leaf);
          if (!text) continue;
          ctx.fillText(text, atX(view.x[leaf]) + 4, atY(view.y[leaf]));
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
      var seen = overview(box);
      var inner = strip.wide - RAIL_PAD * 2;
      var reach = seen.reach;
      var acrossX = (reach.highX - reach.lowX) || 1;
      var atX = function (value) {
        return strip.x0 + RAIL_PAD + ((value - reach.lowX) / acrossX) * inner;
      };
      var atY = function (row) { return heightOfRow(row, box); };

      // The edge it stands behind, so the rail reads as a margin and not as
      // more tree.
      ctx.strokeStyle = theme.frame;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(strip.x0 + 0.5, 0);
      ctx.lineTo(strip.x0 + 0.5, box.tall);
      ctx.stroke();

      ctx.strokeStyle = theme.faint;
      ctx.beginPath();
      for (var each = 0; each < seen.count; each++) {
        var node = seen.nodes[each];
        var up = view.parent[node];
        if (up === 0xffffffff) continue;
        var y = atY(view.y[node]);
        var x1 = atX(view.x[node]);
        var x0 = atX(view.x[up]);
        ctx.moveTo(x0, y);
        ctx.lineTo(x1, y);
        ctx.moveTo(x0, y);
        ctx.lineTo(x0, atY(view.y[up]));
      }
      ctx.stroke();

      // The rows on screen. Three rows out of two million is a thousandth of a
      // pixel, so the mark is held to a size a reader can see and a hand can
      // catch, and it is kept on the rail rather than allowed to hang off it.
      var top = atY(camera.y0);
      var foot = atY(camera.y1);
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
      var span = camera.y1 - camera.y0;
      settle(camera.y0 - (dy / box.tall) * span, span);
    }

    // Is this point, in canvas pixels, on the rail rather than on the tree?
    function onRail(px) {
      var strip = view ? rail(size()) : null;
      return !!strip && px >= strip.x0;
    }

    // Put the rows on screen where this height on the rail points, keeping how
    // many of them there are. A click and a drag are the same gesture: the
    // window follows the hand rather than being nudged by it.
    function scrubTo(py) {
      if (!view) return;
      var box = size();
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
      var from = firstRow(camera.y0);
      var to = firstRow(camera.y1);
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
      onRail: onRail,
      scrubTo: scrubTo,
      goTo: goTo,
      looking: looking,
      loaded: function () { return !!view; },
      count: function () { return view ? view.count : 0; },
      // What the last paint put on the canvas, and the depth it fitted them
      // into. Both are answers about the picture rather than about the tree,
      // which is what a check of the picture needs to ask.
      shown: function () { return view ? view.shown.subarray(0, view.picked || 0) : new Uint32Array(0); },
      depth: function () { return { x0: camera.x0, x1: camera.x1 }; },
    };
  }

  return { make: make };
})();
