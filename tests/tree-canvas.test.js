// The page's canvas, checked by running it. Plain node, no framework and no
// install: the crate has no runtime dependencies and its tests should not
// either. Run it with `node tests/tree-canvas.test.js`.

const fs = require("fs");
const path = require("path");
const assert = require("assert");

// --------------------------------------------------------------- a canvas

// Enough of one for the painter to draw into, and it keeps what it was told to
// draw so a test can ask what is on it.
function fakeCanvas(wide, tall) {
  const strokes = [];
  const rects = [];
  const ctx = {
    setTransform() {}, clearRect() {}, beginPath() {}, stroke() {},
    moveTo(x, y) { strokes.push(["move", x, y, ctx.strokeStyle]); },
    lineTo(x, y) { strokes.push(["line", x, y, ctx.strokeStyle]); },
    fillRect(x, y, w, h) { rects.push({ kind: "fill", x, y, w, h, paint: ctx.fillStyle }); },
    strokeRect(x, y, w, h) { rects.push({ kind: "stroke", x, y, w, h, paint: ctx.strokeStyle }); },
    fillText() {},
    strokeStyle: "", fillStyle: "", lineWidth: 1, font: "", textBaseline: "",
  };
  return {
    width: wide, height: tall, clientWidth: wide, clientHeight: tall,
    getContext: () => ctx,
    strokes,
    rects,
  };
}

global.window = { devicePixelRatio: 1 };
require(path.join(__dirname, "..", "docs", "assets", "tree-canvas.js"));
const canvasModule = global.window.karyonCanvas;

// ----------------------------------------------------------------- a tree

// A balanced tree of `1 << levels` tips, laid out the way `Tree::layout` does:
// a leaf sits on its own row and a parent on the mean of its children's.
function balanced(levels) {
  const tips = 1 << levels;
  const count = 2 * tips - 1;
  const x = new Float32Array(count);
  const y = new Float32Array(count);
  const parent = new Uint32Array(count);
  parent[0] = 0xffffffff;
  const kids = [];
  for (let i = 0; i < count; i++) kids.push([]);
  // Node 0 is the root; each internal node `i` has children `2i+1`, `2i+2`.
  for (let i = 1; i < count; i++) {
    parent[i] = (i - 1) >> 1;
    kids[parent[i]].push(i);
  }
  for (let i = 0; i < count; i++) x[i] = Math.floor(Math.log2(i + 1)) * 0.05;
  let row = 0;
  for (let i = count - 1; i >= 0; i--) {
    if (kids[i].length === 0) y[i] = row++;
  }
  // Leaves came out back to front above, so put them in walk order and then
  // fold the internal rows up from them.
  row = 0;
  const stack = [0];
  const order = [];
  while (stack.length) {
    const at = stack.pop();
    order.push(at);
    for (let k = kids[at].length - 1; k >= 0; k--) stack.push(kids[at][k]);
  }
  for (const at of order) if (kids[at].length === 0) y[at] = row++;
  for (let i = count - 1; i >= 0; i--) {
    if (kids[i].length) y[i] = (y[kids[i][0]] + y[kids[i][kids[i].length - 1]]) / 2;
  }
  return {
    count, x, y, parent,
    start: new Uint32Array(count),
    length: new Uint32Array(count),
    names: new Uint8Array(0),
  };
}

const theme = {
  branch: "#000", muted: "#666", font: "sans-serif",
  frame: "#ccc", faint: "#999", window: "rgba(0,0,255,0.15)", edge: "#00f",
};

// Wide enough for the rail, and narrow enough to be refused one.
const WIDE = 900;
const NARROW = 380;

// ---------------------------------------------------------------- the test

// The property: whatever is drawn is a tree. Every branch on screen but the
// root's has the branch it hangs from on screen too, so the picture has a trunk
// instead of being a hedge of loose strokes.
function drawnSet(placed, tall, wide) {
  const canvas = fakeCanvas(wide || WIDE, tall);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  const report = painter.paint(theme);
  return { painter, report, canvas };
}

// A painter that keeps every row it is given, so the test can see what the
// selection looked like before the ancestors were walked back in. This is the
// same stride the real one uses.
function withoutAncestors(placed, tall) {
  const byRow = Uint32Array.from(
    Array.from({ length: placed.count }, (_, i) => i).sort((a, b) => placed.y[a] - placed.y[b])
  );
  const stride = Math.max(1, Math.ceil(placed.count / tall));
  const kept = new Set();
  for (let at = 0; at < placed.count; at += stride) kept.add(byRow[at]);
  let joined = 0;
  for (const node of kept) {
    if (placed.parent[node] !== 0xffffffff && kept.has(placed.parent[node])) joined += 1;
  }
  return { kept: kept.size, joined };
}

let failures = 0;
function check(what, run) {
  try {
    run();
    console.log("ok   " + what);
  } catch (cause) {
    failures += 1;
    console.log("FAIL " + what + "\n     " + cause.message);
  }
}

check("what is drawn hangs together, at every zoom", () => {
  const placed = balanced(14); // 16,384 tips, far taller than any canvas
  for (const tall of [200, 800, 801, 1600]) {
    const { painter, report } = drawnSet(placed, tall);
    assert.ok(report.drawn > 0, `nothing drawn at ${tall} px`);
    // Ask the painter what it put on the canvas by walking the same selection:
    // every stroke pair is a branch, and the branch's own parent must be there.
    const shown = painter.shown();
    const on = new Set(shown);
    let loose = 0;
    for (const node of shown) {
      if (placed.parent[node] === 0xffffffff) continue;
      if (!on.has(placed.parent[node])) loose += 1;
    }
    assert.strictEqual(loose, 0, `${loose} of ${shown.length} branches hang from nothing at ${tall} px`);
    assert.ok(on.has(0), `the root is not drawn at ${tall} px`);
  }
});

check("and the check bites: stepping over rows alone leaves nothing joined", () => {
  const placed = balanced(14);
  const naive = withoutAncestors(placed, 800);
  assert.strictEqual(
    naive.joined, 0,
    `expected the bare stride to draw ${naive.kept} unconnected branches, got ${naive.joined} joined`
  );
});

check("a drag cannot push the tree off the canvas", () => {
  const placed = balanced(12);
  const { painter } = drawnSet(placed, 800);
  for (let i = 0; i < 500; i++) painter.panBy(0, -80);
  assert.ok(painter.paint(theme).drawn > 0, "dragged down until the canvas was empty");
  for (let i = 0; i < 1000; i++) painter.panBy(0, 80);
  assert.ok(painter.paint(theme).drawn > 0, "dragged up until the canvas was empty");
});

check("the depth window holds the root and the deepest tip on screen", () => {
  const placed = balanced(14);
  const { painter } = drawnSet(placed, 800);
  const shown = painter.shown();
  let low = Infinity, high = -Infinity;
  for (const node of shown) {
    if (placed.x[node] < low) low = placed.x[node];
    if (placed.x[node] > high) high = placed.x[node];
  }
  const window = painter.depth();
  assert.ok(window.x0 <= low + 1e-6, "the root is off the left edge");
  assert.ok(window.x1 >= high - 1e-6, "the deepest tip on screen is off the right edge");
});

// ----------------------------------------------------------------- the rail

check("the rail shows the whole tree, whatever the window holds", () => {
  const placed = balanced(14);
  const { painter, report } = drawnSet(placed, 800);
  assert.ok(report.rail, "no rail on a canvas with room for one");
  const wide = report.rail.drawn;
  // Fly all the way in. The rail is the whole tree and must not follow.
  for (let i = 0; i < 60; i++) painter.zoomAt(200, 400, 1.3);
  const close = painter.paint(theme);
  assert.ok(close.drawn < 100, `expected to be deep in, ${close.drawn} branches on screen`);
  assert.strictEqual(close.rail.drawn, wide, "the rail changed with the zoom");
});

check("and the check bites: a rail built from the window would follow it", () => {
  const placed = balanced(14);
  const { painter, report } = drawnSet(placed, 800);
  for (let i = 0; i < 60; i++) painter.zoomAt(200, 400, 1.3);
  const close = painter.paint(theme);
  // What the main view drew is what a window-built rail would have held, and
  // it is nothing like the whole tree, which is the point of the check above.
  assert.notStrictEqual(close.drawn, report.rail.drawn);
});

check("the mark stays big enough to see and to catch", () => {
  const placed = balanced(16); // 65,536 tips
  const { painter } = drawnSet(placed, 800);
  for (let i = 0; i < 80; i++) painter.zoomAt(200, 400, 1.3);
  const report = painter.paint(theme);
  const rows = painter.looking().rows;
  assert.ok(rows < 20, `expected a handful of rows, got ${rows}`);
  assert.ok(
    report.rail.deep >= 4,
    `the mark for ${rows} rows came out ${report.rail.deep.toFixed(3)} px deep`
  );
  assert.ok(report.rail.top >= 0, "the mark hangs off the top of the rail");
  assert.ok(
    report.rail.top + report.rail.deep <= 800 + 1e-6,
    "the mark hangs off the bottom of the rail"
  );
});

check("the mark is drawn, in the colours it was given", () => {
  const placed = balanced(12);
  const { canvas } = drawnSet(placed, 800);
  const filled = canvas.rects.filter((r) => r.kind === "fill");
  assert.strictEqual(filled.length, 1, "expected one filled mark");
  assert.strictEqual(filled[0].paint, theme.window, "the mark is not the window colour");
  assert.ok(canvas.rects.some((r) => r.kind === "stroke" && r.paint === theme.edge), "the mark has no edge");
});

check("a click on the rail puts those rows on screen", () => {
  const placed = balanced(16);
  const { painter } = drawnSet(placed, 800);
  // Zoomed in, so the mark is small enough to be placed rather than clamped.
  for (let i = 0; i < 20; i++) painter.zoomAt(200, 400, 1.3);
  const before = painter.looking().rows;
  for (const py of [120, 400, 600, 750]) {
    painter.scrubTo(py);
    const report = painter.paint(theme);
    const middle = report.rail.top + report.rail.deep / 2;
    assert.ok(
      Math.abs(middle - py) < 3,
      `clicked at ${py} px and the mark came out centred at ${middle.toFixed(1)}`
    );
    assert.ok(
      Math.abs(painter.looking().rows - before) <= 2,
      "the click changed how many rows are shown"
    );
  }
});

check("and the check bites: the rail's two directions are inverses", () => {
  const placed = balanced(16);
  const { painter } = drawnSet(placed, 800);
  for (let i = 0; i < 20; i++) painter.zoomAt(200, 400, 1.3);
  // A scrub to the very top and the very bottom must not land in the same
  // place, which is what a dropped or constant mapping would do.
  painter.scrubTo(0);
  const top = painter.paint(theme).rail.top;
  painter.scrubTo(800);
  const foot = painter.paint(theme).rail.top;
  assert.ok(foot - top > 700, `the whole rail moved the mark only ${(foot - top).toFixed(1)} px`);
});

check("the rail knows what belongs to it", () => {
  const placed = balanced(12);
  const { painter, report } = drawnSet(placed, 800);
  assert.ok(painter.onRail(report.rail.x0 + 2), "a point on the rail was not claimed");
  assert.ok(painter.onRail(WIDE - 1), "the far edge was not claimed");
  assert.ok(!painter.onRail(report.rail.x0 - 2), "a point on the tree was claimed by the rail");
  assert.ok(!painter.onRail(10), "the root end was claimed by the rail");
});

check("a narrow canvas gets no rail, and all of its width", () => {
  const placed = balanced(12);
  const roomy = drawnSet(placed, 800, WIDE);
  const tight = drawnSet(placed, 800, NARROW);
  assert.ok(!tight.report.rail, "a phone width canvas was given a rail");
  assert.ok(!tight.painter.onRail(NARROW - 1), "the rail claims points on a canvas that has none");
  assert.strictEqual(tight.canvas.rects.length, 0, "something was drawn where the rail would be");
  assert.ok(roomy.report.rail, "a wide canvas was refused a rail");
});

process.exit(failures ? 1 : 0);
