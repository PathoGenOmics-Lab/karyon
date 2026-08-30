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

    function paint(theme) {
      if (!view) return { drawn: 0, skipped: 0 };
      var box = size();
      var ctx = canvas.getContext("2d");
      ctx.setTransform(box.ratio, 0, 0, box.ratio, 0, 0);
      ctx.clearRect(0, 0, box.wide, box.tall);

      var spanX = camera.x1 - camera.x0 || 1;
      var spanY = camera.y1 - camera.y0 || 1;
      var sx = box.wide / spanX;
      var sy = box.tall / spanY;
      var atX = function (value) { return (value - camera.x0) * sx; };
      var atY = function (value) { return (value - camera.y0) * sy; };

      var from = firstRow(camera.y0 - 1);
      var to = firstRow(camera.y1 + 1);
      var wanted = to - from;
      // One row per pixel is all a screen can show. Past that they land on each
      // other, so they are stepped over.
      var stride = Math.max(1, Math.ceil(wanted / Math.max(1, box.tall)));

      // The depth axis follows the rows in view rather than being zoomed
      // alongside them. Zooming both narrowed the window in x as well until it
      // held no branches at all: past fifty wheel notches on a million tip tree
      // the canvas went blank while the rows were still there. A tree is read
      // by scrolling through its tips with the whole depth of what is on screen
      // in front of you, which is also what --focus gives a figure.
      var lowX = Infinity, highX = -Infinity;
      for (var scan = from; scan < to; scan += stride) {
        var it = view.byRow[scan];
        if (view.x[it] < lowX) lowX = view.x[it];
        if (view.x[it] > highX) highX = view.x[it];
        var over = view.parent[it];
        if (over !== 0xffffffff && view.x[over] < lowX) lowX = view.x[over];
      }
      if (!(highX > lowX)) {
        var b = view.bounds;
        lowX = b.lowX;
        highX = b.highX;
      }
      var margin = (highX - lowX) * 0.28;
      camera.x0 = lowX - margin * 0.05;
      camera.x1 = highX + margin;
      spanX = camera.x1 - camera.x0 || 1;
      sx = box.wide / spanX;

      ctx.strokeStyle = theme.branch;
      ctx.lineWidth = 1;
      ctx.beginPath();
      var drawn = 0;
      for (var at = from; at < to; at += stride) {
        var node = view.byRow[at];
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
      return { drawn: drawn, skipped: wanted - drawn, labels: labels, stride: stride, rowsInView: wanted };
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

    function panBy(dx, dy) {
      var box = size();
      var byY = (dy / box.tall) * (camera.y1 - camera.y0);
      camera.y0 -= byY;
      camera.y1 -= byY;
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
      goTo: goTo,
      looking: looking,
      loaded: function () { return !!view; },
      count: function () { return view ? view.count : 0; },
    };
  }

  return { make: make };
})();
