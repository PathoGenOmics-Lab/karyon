// The bridge every page runs the program through, checked where it rewrites
// a command. Plain node, no framework and no install, like the canvas beside
// it. Run it with `node tests/karyon-wasm.test.js`.

const path = require("path");
const assert = require("assert");

// The bridge sets itself on `self`, as a page and a worker both have one, and
// asks nothing else of where it runs until a figure is drawn.
globalThis.self = globalThis;
require(path.join(__dirname, "..", "docs", "assets", "karyon-wasm.js"));
const karyon = self.karyon;

// A drag rewrites the place where it stands. It rewrote the first stretch of
// the command spelled like it, so a --shade written before the place, on the
// same stretch, was what moved, and the figure stayed where it was.
assert.strictEqual(
  karyon.retarget("--shade chr1:1,001-2,000=x chr1:1,001-2,000 d.bg", 1501, 2500),
  "--shade chr1:1,001-2,000=x chr1:1,501-2,500 d.bg"
);
assert.strictEqual(
  karyon.panned("--shade chr1:1-1,000 chr1:1-1,000 d.bg", -0.5),
  "--shade chr1:1-1,000 chr1:501-1,500 d.bg"
);

// The rest of the command is left as it was written, quotes and all, and a
// quoted place is rewritten without its quotes.
assert.strictEqual(
  karyon.retarget("chr1:1-1,000 d.bg --title 'a figure'", 101, 1100),
  "chr1:101-1,100 d.bg --title 'a figure'"
);
assert.strictEqual(karyon.retarget("'chr1:1-1,000' d.bg", 101, 1100), "chr1:101-1,100 d.bg");

// A value that looks like a place is never the place.
assert.strictEqual(
  karyon.retarget("--label chr1:5-9 chr1:1-1,000 d.bg", 101, 1100),
  "--label chr1:5-9 chr1:101-1,100 d.bg"
);

console.log("karyon-wasm: the place is rewritten where it stands");
