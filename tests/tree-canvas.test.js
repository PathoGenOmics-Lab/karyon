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
  const texts = [];
  const arcs = [];
  const ctx = {
    setTransform() {}, clearRect() {}, beginPath() {}, stroke() {},
    moveTo(x, y) { strokes.push(["move", x, y, ctx.strokeStyle]); },
    lineTo(x, y) { strokes.push(["line", x, y, ctx.strokeStyle]); },
    fillRect(x, y, w, h) { rects.push({ kind: "fill", x, y, w, h, paint: ctx.fillStyle }); },
    arc(x, y, r, a0, a1) { arcs.push({ x, y, r, a0, a1 }); },
    drawImage() {},
    strokeRect(x, y, w, h) { rects.push({ kind: "stroke", x, y, w, h, paint: ctx.strokeStyle }); },
    fillText(t, x, y) { texts.push({ t, x, y }); },
    measureText(t) { return { width: t.length * 6 }; },
    strokeStyle: "", fillStyle: "", lineWidth: 1, font: "", textBaseline: "",
  };
  return {
    width: wide, height: tall, clientWidth: wide, clientHeight: tall,
    getContext: () => ctx,
    strokes,
    rects,
    texts,
    arcs,
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
  plate: "#fff",
};

// Wide enough for the rail, narrow enough to get only a thin one, and too
// narrow to be given one at all.
const WIDE = 900;
const SNUG = 380;
const NARROW = 220;

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
    // Every branch on screen hangs from a branch on screen, or from the root.
    // A parent stepped over by the pixel test is not drawn as a branch of its
    // own, so what has to be there is the end the branch was redirected to.
    const shown = painter.shown();
    const ends = new Set();
    for (let i = 0; i < shown.count; i++) ends.add(shown.child[i]);
    let loose = 0;
    for (let i = 0; i < shown.count; i++) {
      const over = shown.up[i];
      if (placed.parent[over] === 0xffffffff) continue;
      if (!ends.has(over)) loose += 1;
    }
    assert.strictEqual(loose, 0, `${loose} of ${shown.count} branches hang from nothing at ${tall} px`);
    let root = false;
    for (let i = 0; i < shown.count; i++) {
      if (placed.parent[shown.up[i]] === 0xffffffff) root = true;
    }
    assert.ok(root, `nothing reaches the root at ${tall} px`);
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
  for (let i = 0; i < shown.count; i++) {
    for (const node of [shown.child[i], shown.up[i]]) {
      if (placed.x[node] < low) low = placed.x[node];
      if (placed.x[node] > high) high = placed.x[node];
    }
  }
  const window = painter.depth();
  assert.ok(window.x0 <= low + 1e-6, "the root is off the left edge");
  assert.ok(window.x1 >= high - 1e-6, "the deepest tip on screen is off the right edge");
});

// A ladder is the shape that used to make the walk to the root cost the whole
// tree on every frame, so it is the shape the bound is checked on.
function ladder(tips) {
  const count = 2 * tips - 1;
  const x = new Float32Array(count);
  const y = new Float32Array(count);
  const parent = new Uint32Array(count);
  parent[0] = 0xffffffff;
  // Node 0 is the root. Each internal node i has a tip and the next internal.
  let spine = 0;
  let row = 0;
  for (let i = 1; i < count; i += 2) {
    const tip = i, next = i + 1;
    parent[tip] = spine;
    x[tip] = x[spine] + 0.01;
    y[tip] = row++;
    if (next < count) {
      parent[next] = spine;
      x[next] = x[spine] + 0.005;
      spine = next;
    }
  }
  // The spine sits on the mean of what hangs below it, which for a ladder is
  // near enough the middle of the rows still to come.
  for (let i = count - 1; i >= 0; i--) if (!y[i] && i) y[i] = row / 2;
  return { count, x, y, parent, start: new Uint32Array(count), length: new Uint32Array(count), names: new Uint8Array(0) };
}

check("a ladder cannot make one frame walk the whole tree", () => {
  const placed = ladder(60000); // 119,999 nodes on one spine
  const { painter, report } = drawnSet(placed, 800);
  assert.ok(
    report.drawn < 20000,
    `a ladder of ${placed.count} nodes drew ${report.drawn} branches in one frame`
  );
  for (let i = 0; i < 40; i++) painter.zoomAt(200, 400, 1.3);
  const close = painter.paint(theme);
  assert.ok(
    close.drawn < 20000,
    `zoomed in on a ladder it still drew ${close.drawn} branches`
  );
});

check("and the check bites: without the pixel step a ladder walks all of it", () => {
  const placed = ladder(60000);
  // What the walk would cost with no step over ancestors that share a pixel:
  // every tip hangs off the spine, so one walk is half the tree.
  const byRow = Uint32Array.from(
    Array.from({ length: placed.count }, (_, i) => i).sort((a, b) => placed.y[a] - placed.y[b])
  );
  const seen = new Uint8Array(placed.count);
  let reached = 0;
  for (let at = 0; at < placed.count; at += Math.ceil(placed.count / 800)) {
    let walk = byRow[at];
    while (walk !== 0xffffffff && !seen[walk]) { seen[walk] = 1; reached += 1; walk = placed.parent[walk]; }
  }
  assert.ok(reached > 50000, `expected the bare walk to reach most of the tree, it reached ${reached}`);
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
    painter.jumpTo(0, py);
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
  painter.jumpTo(0, 0);
  const top = painter.paint(theme).rail.top;
  painter.jumpTo(0, 800);
  const foot = painter.paint(theme).rail.top;
  assert.ok(foot - top > 700, `the whole rail moved the mark only ${(foot - top).toFixed(1)} px`);
});

check("the rail knows what belongs to it", () => {
  const placed = balanced(12);
  const { painter, report } = drawnSet(placed, 800);
  assert.ok(painter.onMap(report.rail.x0 + 2, 400), "a point on the rail was not claimed");
  assert.ok(painter.onMap(WIDE - 1, 400), "the far edge was not claimed");
  assert.ok(!painter.onMap(report.rail.x0 - 2, 400), "a point on the tree was claimed by the rail");
  assert.ok(!painter.onMap(10, 400), "the root end was claimed by the rail");
});

check("a narrow canvas gets a thin rail, and a tiny one gets none", () => {
  const placed = balanced(12);
  const roomy = drawnSet(placed, 800, WIDE);
  const snug = drawnSet(placed, 800, SNUG);
  const tight = drawnSet(placed, 800, NARROW);

  assert.ok(roomy.report.rail, "a wide canvas was refused a rail");
  assert.ok(snug.report.rail, "a phone width canvas was refused a rail");
  assert.ok(
    snug.report.rail.wide < roomy.report.rail.wide,
    `the rail did not narrow with the canvas: ${snug.report.rail.wide} against ${roomy.report.rail.wide}`
  );
  assert.ok(snug.report.rail.wide >= 30, "the thin rail is too thin to read or to catch");

  assert.ok(!tight.report.rail, "a canvas with no room at all was given a rail");
  assert.ok(!tight.painter.onMap(NARROW - 1, 400), "the rail claims points on a canvas that has none");
  assert.strictEqual(tight.canvas.rects.length, 0, "something was drawn where the rail would be");
});

check("and the check bites: a rail of a fixed width would not narrow", () => {
  const placed = balanced(12);
  const roomy = drawnSet(placed, 800, WIDE).report.rail;
  const snug = drawnSet(placed, 800, SNUG).report.rail;
  // The wide canvas is held at the cap, so if the rule were a constant these
  // two would be equal and the check above would be about nothing.
  assert.notStrictEqual(snug.wide, roomy.wide);
  assert.ok(roomy.wide <= 78, "the rail grew past what it is capped at");
});

// ----------------------------------------------------------------- the disc

// The crate's own arithmetic, written out here from radial.rs rather than
// borrowed from the module under test, so the check is a second opinion and
// not an echo. A row of `terminals` goes to an angle starting at the top of
// the circle and going all the way round; a depth goes to a fraction of the
// way from a hole eight percent out to the rim.
function crateAngle(row, lowY, terminals) {
  return -Math.PI / 2 + (Math.PI * 2 * (row - lowY)) / terminals;
}
function crateRadius(depth, lowX, highX) {
  const span = highX - lowX;
  const part = span > 0 ? (depth - lowX) / span : 0;
  return 0.08 + Math.min(1, Math.max(0, part)) * 0.92;
}

check("the disc puts a node where the crate would", () => {
  const placed = balanced(10); // 1,024 tips
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.shape("disc");
  painter.paint(theme);

  let lowY = Infinity, highY = -Infinity, lowX = Infinity, highX = -Infinity;
  for (let i = 0; i < placed.count; i++) {
    if (placed.y[i] < lowY) lowY = placed.y[i];
    if (placed.y[i] > highY) highY = placed.y[i];
    if (placed.x[i] < lowX) lowX = placed.x[i];
    if (placed.x[i] > highX) highX = placed.x[i];
  }
  const terminals = highY - lowY + 1;
  // At Fit the camera is the whole disc: half a canvas height over 1.1.
  const scale = 800 / 2.2;
  const midX = WIDE / 2, midY = 800 / 2;

  let worst = 0;
  for (const node of [0, 1, 2, 17, 500, placed.count - 1]) {
    const angle = crateAngle(placed.y[node], lowY, terminals);
    const radius = crateRadius(placed.x[node], lowX, highX) * scale;
    const want = { x: midX + Math.cos(angle) * radius, y: midY + Math.sin(angle) * radius };
    const got = painter.where(node);
    worst = Math.max(worst, Math.abs(got.x - want.x), Math.abs(got.y - want.y));
  }
  assert.ok(worst < 0.001, `the canvas and the crate disagree by ${worst.toFixed(4)} px`);
});

check("and the check bites: the disc is not the rectangle in disguise", () => {
  const placed = balanced(10);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  const flat = painter.where(3);
  painter.shape("disc");
  painter.paint(theme);
  const bent = painter.where(3);
  assert.ok(
    Math.abs(flat.x - bent.x) > 1 || Math.abs(flat.y - bent.y) > 1,
    "switching projection moved nothing"
  );
});

check("switching projection asks the program for nothing", () => {
  const placed = balanced(10);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  const before = { x: placed.x.slice(), y: placed.y.slice(), parent: placed.parent.slice() };
  painter.shape("disc");
  painter.paint(theme);
  painter.shape("rows");
  painter.paint(theme);
  // The layout is the same rows and depths whichever way it is drawn, so
  // nothing here may have touched them.
  assert.deepStrictEqual(Array.from(placed.x), Array.from(before.x), "the depths moved");
  assert.deepStrictEqual(Array.from(placed.y), Array.from(before.y), "the rows moved");
  assert.deepStrictEqual(Array.from(placed.parent), Array.from(before.parent), "the tree moved");
});

check("the walk gives up at the edge of the canvas", () => {
  const placed = balanced(15); // 32,768 tips
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.shape("disc");
  const whole = painter.paint(theme).drawn;
  // In on the rim, where the middle of the disc that every walk heads for is
  // off the canvas.
  for (let i = 0; i < 12; i++) painter.zoomAt(WIDE / 2 + 200, 400, 1.35);
  const rim = painter.paint(theme).drawn;
  assert.ok(rim < whole, `zoomed in on the rim it drew ${rim}, more than the ${whole} of the whole disc`);
  assert.ok(rim > 0, "zoomed in on the rim it drew nothing at all");
});

check("and the check bites: without the edge the walk reaches the middle", () => {
  const placed = balanced(15);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.shape("disc");
  painter.paint(theme);
  for (let i = 0; i < 12; i++) painter.zoomAt(WIDE / 2 + 200, 400, 1.35);
  painter.paint(theme);
  // The root is what every walk is heading for, and at this zoom it is far
  // outside the canvas: if the walk did not stop, it would be drawn.
  const seen = painter.shown();
  let root = 0;
  for (let i = 0; i < placed.count; i++) if (placed.parent[i] === 0xffffffff) root = i;
  const at = painter.where(root);
  assert.ok(
    at.x < -48 || at.x > WIDE + 48 || at.y < -48 || at.y > 848,
    "the root is on the canvas, so this check proves nothing"
  );
  let drawnRoot = false;
  for (let i = 0; i < seen.count; i++) if (seen.child[i] === root) drawnRoot = true;
  assert.ok(!drawnRoot, "the root was drawn although it is off the canvas");
});

check("the dial is the whole disc, and it does not follow the zoom", () => {
  const placed = balanced(12);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.shape("disc");
  const wide = painter.paint(theme).rail;
  assert.ok(wide, "no dial on a canvas with room for one");
  for (let i = 0; i < 20; i++) painter.zoomAt(WIDE / 2, 400, 1.3);
  const close = painter.paint(theme).rail;
  assert.strictEqual(close.drawn, wide.drawn, "the dial changed with the zoom");
  assert.ok(close.wide <= wide.wide, "the window on the dial grew as the view shrank");
  assert.ok(close.deep >= 4 && close.wide >= 4, "the window on the dial is too small to catch");
});

check("a click on the dial goes to that part of the disc", () => {
  const placed = balanced(12);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.shape("disc");
  painter.paint(theme);
  for (let i = 0; i < 14; i++) painter.zoomAt(WIDE / 2, 400, 1.3);
  const spot = painter.paint(theme).rail;
  const midX = spot.x0 + spot.side / 2, midY = spot.y0 + spot.side / 2;
  for (const [dx, dy] of [[-30, -30], [30, -30], [30, 30], [-30, 30]]) {
    painter.jumpTo(midX + dx, midY + dy);
    const after = painter.paint(theme).rail;
    const seen = [after.left + after.wide / 2 - midX, after.top + after.deep / 2 - midY];
    assert.ok(
      Math.abs(seen[0] - dx) < 4 && Math.abs(seen[1] - dy) < 4,
      `clicked ${dx},${dy} from the middle and the window went to ${seen[0].toFixed(1)},${seen[1].toFixed(1)}`
    );
  }
});

check("the dial knows what belongs to it, and a narrow canvas gets none", () => {
  const placed = balanced(12);
  const roomy = canvasModule.make(fakeCanvas(WIDE, 800));
  roomy.load(placed);
  roomy.shape("disc");
  const spot = roomy.paint(theme).rail;
  assert.ok(roomy.onMap(spot.x0 + 4, spot.y0 + 4), "a point on the dial was not claimed");
  assert.ok(!roomy.onMap(spot.x0 - 4, spot.y0 - 4), "a point on the disc was claimed by the dial");
  assert.ok(!roomy.onMap(WIDE / 2, 400), "the middle of the disc was claimed by the dial");

  const snug = canvasModule.make(fakeCanvas(SNUG, 800));
  snug.load(placed);
  snug.shape("disc");
  const small = snug.paint(theme).rail;
  assert.ok(small, "a phone width canvas was refused a dial");
  assert.ok(small.side < spot.side, "the dial did not shrink with the canvas");
  assert.ok(small.side <= SNUG * 0.34 + 1, "the dial takes more than a third of the width");

  const tight = canvasModule.make(fakeCanvas(NARROW, 800));
  tight.load(placed);
  tight.shape("disc");
  assert.ok(!tight.paint(theme).rail, "a canvas with no room at all was given a dial");
  assert.ok(!tight.onMap(NARROW - 4, 796), "the dial claims points on a canvas that has none");
});

// --------------------------------------------------------------- rootless

// A layout with no root, in the shape the wire delivers one: positions in a
// plane, what each branch hangs from once the tree is re-rooted at its middle,
// and the order the terminals come round it. Two rings of tips off a middle,
// so it has structure rather than being a wheel of spokes.
function spread(spokes) {
  const count = 1 + spokes * 2;
  const x = new Float32Array(count);
  const y = new Float32Array(count);
  const parent = new Uint32Array(count);
  parent[0] = 0xffffffff;
  const order = new Uint32Array(spokes);
  for (let i = 0; i < spokes; i++) {
    const inner = 1 + i * 2;
    const tip = inner + 1;
    const turn = (Math.PI * 2 * i) / spokes;
    x[inner] = Math.cos(turn) * (0.4 + (i % 3) * 0.05);
    y[inner] = Math.sin(turn) * (0.4 + (i % 3) * 0.05);
    x[tip] = Math.cos(turn) * (1 + (i % 5) * 0.1);
    y[tip] = Math.sin(turn) * (1 + (i % 5) * 0.1);
    parent[inner] = 0;
    parent[tip] = inner;
    order[i] = tip;
  }
  return {
    count, x, y, parent, order,
    start: new Uint32Array(count),
    length: new Uint32Array(count),
    names: new Uint8Array(0),
  };
}

function rootlessPainter(spokes, wide, tall) {
  const placed = spread(spokes || 600);
  const canvas = fakeCanvas(wide || WIDE, tall || 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  return { painter, placed, canvas };
}

// ------------------------------------------------------------- the search

// A balanced tree whose tips carry names, since the search is about names and
// the fixture above has none.
function named(levels) {
  const placed = balanced(levels);
  const words = [];
  const start = new Uint32Array(placed.count);
  const length = new Uint32Array(placed.count);
  let blob = "";
  for (let node = 0; node < placed.count; node++) {
    // A leaf is a node nothing points to as a parent, and in this fixture the
    // back half of the array is exactly that.
    if (node < (placed.count - 1) / 2) continue;
    const lineage = node % 3 === 0 ? "L2" : node % 3 === 1 ? "L4.9" : "L1";
    const text = lineage + "_" + String(node).padStart(5, "0");
    start[node] = blob.length;
    length[node] = text.length;
    blob += text;
    words.push(text);
  }
  placed.start = start;
  placed.length = length;
  placed.names = new TextEncoder().encode(blob);
  return { placed, words };
}

function searcher(levels) {
  const { placed, words } = named(levels || 10);
  const painter = canvasModule.make(fakeCanvas(WIDE, 800));
  painter.load(placed);
  painter.paint(theme);
  return { painter, placed, words };
}

check("the search answers with every tip that matches, not the first", () => {
  const { painter, words } = searcher(10);
  const lineage = words.filter((w) => w.indexOf("L4.9") === 0).length;
  assert.ok(lineage > 10, `the fixture has only ${lineage} tips in that lineage`);
  const hits = painter.find("L4.9", "starts");
  assert.strictEqual(hits.length, lineage, `found ${hits.length} of ${lineage}`);
  // And one name is one hit.
  assert.strictEqual(painter.find(words[0], "exact").length, 1);
});

check("and the check bites: a search that stopped at the first would find one", () => {
  const { painter } = searcher(10);
  assert.ok(painter.find("L4.9", "starts").length > 1, "the set has one thing in it");
});

check("it ignores case, takes a part of a name, and takes a list", () => {
  const { painter, words } = searcher(10);
  const one = words[0];
  assert.strictEqual(painter.find(one.toLowerCase(), "loose").length, 1, "case was not folded");
  assert.strictEqual(painter.find(one, "exact").length, 1);
  // The number in the middle of a name, which no prefix search would find.
  const middle = one.slice(3);
  assert.ok(painter.find(middle, "in").length >= 1, "a part of a name found nothing");
  assert.strictEqual(painter.find(middle, "starts").length, 0, "a prefix search matched the middle");
  // A list, separated however it was pasted.
  const list = painter.find(words[0] + ", " + words[1] + "\n" + words[2], "exact");
  assert.strictEqual(list.length, 3, `a list of three found ${list.length}`);
});

check("what it finds comes back in row order", () => {
  const { painter, placed } = searcher(10);
  const hits = painter.find("L", "starts");
  assert.ok(hits.length > 20, "not enough hits to say anything about their order");
  for (let i = 1; i < hits.length; i++) {
    assert.ok(
      placed.y[hits[i]] >= placed.y[hits[i - 1]],
      `hit ${i} sits above the one before it`
    );
  }
});

check("Enter again moves to the next of them, and comes round", () => {
  const { painter, words } = searcher(10);
  assert.ok(painter.goTo("L4.9", "starts"), "the search found nothing");
  const many = painter.found().count;
  assert.ok(many > 3, `only ${many} to step through`);
  assert.strictEqual(painter.found().at, 0);
  assert.deepStrictEqual(painter.nextFound(1), { at: 1, of: many });
  assert.deepStrictEqual(painter.nextFound(1), { at: 2, of: many });
  assert.deepStrictEqual(painter.nextFound(-1), { at: 1, of: many });
  // Round the end from the start.
  painter.nextFound(-1);
  assert.strictEqual(painter.found().at, 0);
  assert.strictEqual(painter.nextFound(-1).at, many - 1, "it did not come round");
  assert.ok(words.length > 0);
});

// ------------------------------------------------------------- the strips

// A layout that arrives with its trait columns already resolved, in the shape
// the wire delivers them: levels with a colour for each scheme, and one level
// index per node.
function striped(levels) {
  // Named tips, because a strip is drawn against a row that has a name: an
  // unnamed internal node is a branch, not a sample.
  const placed = named(10).placed;
  const many = levels || 3;
  const of = new Uint32Array(placed.count);
  for (let node = 0; node < placed.count; node++) {
    // A quarter of the nodes carry nothing, which is what a sheet that names
    // only some of the tips looks like.
    of[node] = node % 4 === 3 ? 0xffffffff : node % many;
  }
  // Two sets that share nothing, so a check can tell which scheme was used.
  const pale = ["#0072b2", "#009e73", "#d55e00", "#7b3294", "#b34f86"];
  const deep = ["#1a1a1a", "#2b2b2b", "#3c3c3c", "#4d4d4d", "#5e5e5e"];
  placed.strips = [
    {
      key: "place",
      label: "place",
      levels: Array.from({ length: many }, (_, at) => ({
        value: "place " + at,
        light: pale[at % pale.length],
        dark: deep[at % deep.length],
      })),
      of,
    },
  ];
  return placed;
}

check("a column of traits is drawn beside the rows, in the colours it arrived with", () => {
  const placed = striped(3);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  const report = painter.paint(theme);
  assert.ok(report.cells > 0, "no trait cells were drawn");

  const cells = canvas.rects.filter((r) => r.kind === "fill");
  const inks = new Set(cells.map((r) => r.paint));
  for (const level of placed.strips[0].levels) {
    assert.ok(inks.has(level.light), `the colour ${level.light} was never used`);
  }
  // And nothing was invented.
  const allowed = new Set(placed.strips[0].levels.map((l) => l.light).concat([theme.window]));
  for (const cell of cells) {
    assert.ok(allowed.has(cell.paint), `a cell was painted ${cell.paint}, which came from nowhere`);
  }
});

check("and the check bites: a tree with no columns draws no cells", () => {
  const placed = balanced(10);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  const report = painter.paint(theme);
  assert.strictEqual(report.cells, 0, "cells were drawn for a tree that has no columns");
});

check("a node the sheet says nothing about gets no cell", () => {
  const placed = striped(3);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  const report = painter.paint(theme);

  // Exactly the named rows that carry a value, and not one more: at this zoom
  // every row is in view, so the count is the whole of the answer.
  let wanted = 0, silent = 0;
  for (let node = 0; node < placed.count; node++) {
    if (!placed.length[node]) continue;
    if (placed.strips[0].of[node] === 0xffffffff) silent += 1;
    else wanted += 1;
  }
  assert.ok(silent > 0, "the fixture has nothing the sheet is silent about");
  assert.strictEqual(report.cells, wanted, `${report.cells} cells for ${wanted} values`);
});

check("the columns take their width from the tree, not from the rail", () => {
  const bareTree = named(10).placed;
  const bare = drawnSet(bareTree, 800, WIDE);
  const withColumns = (() => {
    const painter = canvasModule.make(fakeCanvas(WIDE, 800));
    painter.load(striped(3));
    painter.paint(theme);
    return painter;
  })();
  // The tree is drawn into what is left once the strips have taken their room,
  // so the deepest branch lands further left than it does without them.
  const deep = (p, placed) => {
    let node = 0;
    for (let i = 0; i < placed.count; i++) if (placed.x[i] > placed.x[node]) node = i;
    return p.where(node).x;
  };
  const plain = deep(bare.painter, bareTree);
  const shifted = deep(withColumns, striped(3));
  assert.ok(shifted < plain, `the strips took no room: ${shifted} against ${plain}`);
});

check("the dark scheme uses the dark colour it was given", () => {
  const placed = striped(2);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.paint(Object.assign({}, theme, { dark: true }));
  const inks = new Set(canvas.rects.filter((r) => r.kind === "fill").map((r) => r.paint));
  for (const level of placed.strips[0].levels) {
    assert.ok(inks.has(level.dark), `the dark colour ${level.dark} was never used`);
    assert.ok(!inks.has(level.light), `the light colour ${level.light} was used in the dark`);
  }
});

// --------------------------------------------------------------- the hand

check("the canvas says what is under a point, and nothing where there is none", () => {
  const placed = balanced(12); // 4,096 tips
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.paint(theme);

  // Ask about a branch by asking where one is first, so the check is about the
  // answer and not about guessing a coordinate.
  const shown = painter.shown();
  const wanted = shown.child[Math.floor(shown.count / 2)];
  const spot = painter.where(wanted);
  const found = painter.at(spot.x, spot.y);
  assert.ok(found, "nothing was found where a branch was drawn");
  assert.strictEqual(found.node, wanted, "it named a different branch");
  assert.ok(found.tips >= 1, `a branch with ${found.tips} tips beyond it`);

  // Far from anything drawn.
  assert.strictEqual(painter.at(-500, -500), null, "it found something off the canvas");
});

check("and the check bites: it does not answer for whatever is nearest", () => {
  const placed = balanced(12);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.paint(theme);
  const shown = painter.shown();
  const wanted = shown.child[10];
  const spot = painter.where(wanted);
  // Half the canvas away from any branch, so an answer here would mean it
  // returns the nearest thing however far off it is.
  assert.strictEqual(painter.at(spot.x + 400, spot.y + 400) === null ||
    painter.at(spot.x + 400, spot.y + 400).node !== wanted, true);
  assert.strictEqual(painter.at(spot.x, spot.y).node, wanted, "and it still finds the near one");
});

check("how many tips lie beyond a branch is what the tree says", () => {
  const placed = balanced(10); // 1,024 tips, so the answers are round
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.paint(theme);
  const root = (() => { for (let i = 0; i < placed.count; i++) if (placed.parent[i] === 0xffffffff) return i; })();
  const spot = painter.where(root);
  const found = painter.at(spot.x, spot.y);
  assert.ok(found, "the root was not found where it was drawn");
  assert.strictEqual(found.tips, 1024, `the root has ${found.tips} tips beyond it, not 1024`);
  // Each of the root's two children carries half of them.
  const half = painter.at(painter.where(1).x, painter.where(1).y);
  assert.ok(half && half.tips === 512, `a child of the root has ${half && half.tips}, not 512`);
});

check("taking a clade puts that clade on the screen", () => {
  const placed = balanced(13); // 8,192 tips
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.paint(theme);
  const whole = painter.looking().rows;
  // Node 3 is a quarter of the tree in a balanced one.
  const found = painter.at(painter.where(3).x, painter.where(3).y) || { node: 3, tips: 0 };
  painter.focusOn(3);
  painter.paint(theme);
  const close = painter.looking().rows;
  assert.ok(close < whole / 3, `taking a quarter of the tree left ${close} of ${whole} rows`);
  assert.ok(close > 1, "it took nothing at all");
});

check("and the check bites: taking the root changes nothing", () => {
  const placed = balanced(13);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.paint(theme);
  const whole = painter.looking().rows;
  let root = 0;
  for (let i = 0; i < placed.count; i++) if (placed.parent[i] === 0xffffffff) root = i;
  painter.focusOn(root);
  painter.paint(theme);
  // The root's clade is the tree, so this is the control: a focus that reported
  // a big change here would be reporting the zoom and not the clade.
  assert.ok(
    Math.abs(painter.looking().rows - whole) < whole * 0.3,
    `taking the whole tree went from ${whole} rows to ${painter.looking().rows}`
  );
});

// -------------------------------------------------------------- the wheel

// One wheel notch as the page makes it, from `Math.exp(-deltaY * 0.002)`.
const NOTCH = Math.exp(0.2);

check("a wheel at the middle of a disc never empties it", () => {
  const placed = balanced(14);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.shape("disc");
  let least = Infinity;
  for (let n = 0; n < 40; n++) {
    painter.zoomAt(WIDE / 2, 400, NOTCH);
    const drawn = painter.paint(theme).drawn;
    if (drawn < least) least = drawn;
    assert.ok(drawn > 0, `the canvas came up empty after ${n + 1} notches`);
  }
  assert.ok(least < 400, `it never got close in: the fewest branches was ${least}`);
});

check("and the check bites: the middle of a disc is a hole to fall into", () => {
  // The crate leaves the middle eight per cent of the radius empty, and the
  // camera starts 1.1 units tall. Turning the wheel in the middle without
  // holding on to anything shrinks a window centred on nothing, and this is
  // the notch at which every corner of it is inside that hole.
  const wide = 1.1 * (WIDE / 800);
  let fell = 0;
  for (let n = 1; n <= 40; n++) {
    const half = 1.1 / Math.pow(NOTCH, n);
    const across = wide / Math.pow(NOTCH, n);
    if (Math.sqrt(half * half + across * across) < 0.08) { fell = n; break; }
  }
  assert.ok(fell > 0 && fell < 40, `a plain zoom never falls into the hole, so the check above is empty`);
});

check("a wheel at the middle of a disc keeps narrowing what is in view", () => {
  const placed = balanced(14);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  painter.shape("disc");
  painter.paint(theme);
  const wide = painter.looking().rows;
  for (let n = 0; n < 16; n++) painter.zoomAt(WIDE / 2, 400, NOTCH);
  painter.paint(theme);
  const close = painter.looking().rows;
  assert.ok(
    close < wide / 2,
    `sixteen notches took ${wide} rows to ${close}, which is not going in`
  );
});

check("a wheel cannot take the window off the drawing", () => {
  const { painter } = rootlessPainter(600);
  painter.paint(theme);
  // At the very edge of the canvas, where the drawing is behind you rather
  // than under you.
  for (let n = 0; n < 30; n++) {
    painter.zoomAt(WIDE - 45, 400, NOTCH);
    painter.paint(theme);
  }
  const where = painter.depth();
  assert.ok(where, "the camera came back as nothing");
  // Something of the tree is still within reach of the window: the check is on
  // the outline, so a gap between two branches is still a gap.
  let seen = 0;
  for (let node = 0; node < 1200; node++) {
    const at = painter.where(node);
    if (at.x > -600 && at.x < WIDE + 600 && at.y > -600 && at.y < 1400) seen += 1;
  }
  assert.ok(seen > 0, "the drawing ended up entirely behind the window");
});

check("a wheel into a gap between branches does not empty the canvas either", () => {
  // The corners of a rootless drawing are the gaps between its spokes, so a
  // window can sit squarely on the drawing and hold none of it. Checking the
  // outline is not enough there; what the last paint drew has to be asked.
  for (const [ax, ay] of [[WIDE * 0.9, 80], [WIDE * 0.1, 720]]) {
    const { painter } = rootlessPainter(600);
    painter.paint(theme);
    for (let n = 0; n < 30; n++) {
      painter.zoomAt(ax, ay, NOTCH);
      // What a reader does: turn, look, and if there is nothing there step
      // back. The page does the same in `repaint`.
      let drawn = painter.paint(theme).drawn;
      if (!drawn && painter.stepBack()) drawn = painter.paint(theme).drawn;
      assert.ok(
        drawn > 0,
        `zooming at ${Math.round(ax)},${ay} emptied the canvas after ${n + 1} notches`
      );
    }
  }
});

check("the rows a rootless view holds are its tips and not its coordinates", () => {
  const { painter, placed } = rootlessPainter(600);
  painter.paint(theme);
  // Reading `bounds` as row numbers there gave 1, and the figure the page then
  // asked the program for was folded to eight rows.
  assert.strictEqual(painter.looking().rows, 600, "the tips were not counted");
  assert.notStrictEqual(painter.looking().rows, 1);
});

check("a rootless layout comes up rootless, and stays that way", () => {
  const { painter } = rootlessPainter();
  assert.ok(painter.rootless(), "the layout was not recognised as rootless");
  assert.strictEqual(painter.shapeNow(), "spread", "it did not open in the projection it is for");
  // The other two read rows and depths, which this layout does not have.
  assert.strictEqual(painter.shape("rows"), "spread", "it agreed to be drawn as rows");
  assert.strictEqual(painter.shape("disc"), "spread", "it agreed to be drawn as a disc");
});

check("and the check bites: a rooted layout refuses the rootless projection", () => {
  const placed = balanced(10);
  const canvas = fakeCanvas(WIDE, 800);
  const painter = canvasModule.make(canvas);
  painter.load(placed);
  assert.ok(!painter.rootless());
  assert.strictEqual(painter.shape("spread"), "rows", "rows agreed to be drawn without a root");
  assert.strictEqual(painter.shape("disc"), "disc", "and the one it can do was refused");
});

check("the drawing is the layout to scale, and not squashed", () => {
  const { painter, placed } = rootlessPainter();
  painter.paint(theme);
  // A projection with no root is a similarity: every distance on the canvas is
  // the distance in the layout times one number, the same number in both
  // directions. Anything else has bent the tree.
  const pairs = [[0, 2], [2, 4], [4, 200], [200, 601], [1, 999], [3, 1199]];
  const ratios = [];
  for (const [a, b] of pairs) {
    const here = painter.where(a), there = painter.where(b);
    const drawn = Math.hypot(here.x - there.x, here.y - there.y);
    const laid = Math.hypot(placed.x[a] - placed.x[b], placed.y[a] - placed.y[b]);
    if (laid < 1e-9) continue;
    ratios.push(drawn / laid);
  }
  assert.ok(ratios.length >= 5, "not enough pairs to say anything");
  const low = Math.min(...ratios), high = Math.max(...ratios);
  assert.ok(
    (high - low) / high < 1e-6,
    `the same distance came out ${low.toFixed(3)} to ${high.toFixed(3)} times its size`
  );
});

check("and the check bites: one axis at a different scale fails it", () => {
  const { painter, placed } = rootlessPainter();
  painter.paint(theme);
  const a = painter.where(2), b = painter.where(2 + Math.floor(600 / 4) * 2);
  // Two tips a quarter of the way round from each other are mostly apart in
  // different directions, so a squashed drawing would show up between them.
  assert.ok(Math.abs(a.x - b.x) > 1 && Math.abs(a.y - b.y) > 1, "these two do not test both axes");
});

check("the dial holds the whole tree, and the window shrinks into it", () => {
  const { painter } = rootlessPainter();
  const wide = painter.paint(theme).rail;
  assert.ok(wide, "no dial on a canvas with room for one");
  for (let i = 0; i < 16; i++) painter.zoomAt(WIDE / 2, 400, 1.3);
  const close = painter.paint(theme).rail;
  assert.strictEqual(close.drawn, wide.drawn, "the dial changed with the zoom");
  assert.ok(close.wide < wide.wide && close.deep < wide.deep, "the window did not shrink");
  assert.ok(close.wide >= 4 && close.deep >= 4, "the window is too small to catch");
});

check("nothing is drawn that could not cross the canvas", () => {
  const { painter, placed } = rootlessPainter();
  painter.paint(theme);
  for (let i = 0; i < 10; i++) painter.zoomAt(WIDE / 2 + 260, 300, 1.3);
  painter.paint(theme);
  const seen = painter.shown();
  const code = (q) => {
    let out = 0;
    if (q.x < -24) out |= 1;
    if (q.x > WIDE + 24) out |= 2;
    if (q.y < -24) out |= 4;
    if (q.y > 824) out |= 8;
    return out;
  };
  let impossible = 0;
  for (let i = 0; i < seen.count; i++) {
    if (code(painter.where(seen.child[i])) & code(painter.where(seen.up[i]))) impossible += 1;
  }
  assert.strictEqual(impossible, 0, `${impossible} of ${seen.count} branches cannot reach the canvas`);
  assert.ok(seen.count > 0, "nothing was drawn at all");
});

check("and the check bites: the whole tree has branches that cannot reach it", () => {
  const { painter } = rootlessPainter();
  painter.paint(theme);
  for (let i = 0; i < 10; i++) painter.zoomAt(WIDE / 2 + 260, 300, 1.3);
  painter.paint(theme);
  const code = (q) => {
    let out = 0;
    if (q.x < -24) out |= 1;
    if (q.x > WIDE + 24) out |= 2;
    if (q.y < -24) out |= 4;
    if (q.y > 824) out |= 8;
    return out;
  };
  // Over the whole tree, plenty of branches are off one side together. If the
  // count above were zero for that reason rather than because they were
  // dropped, this would be zero too.
  let offscreen = 0;
  const placed = spread(600);
  for (let node = 1; node < placed.count; node++) {
    const over = placed.parent[node];
    if (over === 0xffffffff) continue;
    if (code(painter.where(node)) & code(painter.where(over))) offscreen += 1;
  }
  assert.ok(offscreen > 100, `only ${offscreen} branches are off one side of the canvas together`);
});

process.exit(failures ? 1 : 0);
