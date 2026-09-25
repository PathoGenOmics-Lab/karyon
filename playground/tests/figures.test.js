// The committed figures, drawn by the wasm a page loads, from outside it. Plain
// node, no framework and no install, like the page's own checks under tests/.
//
// Build the module first, the way CI does, and then run this:
//
//   cargo build --release --target wasm32-unknown-unknown \
//     --manifest-path playground/Cargo.toml
//   node playground/tests/figures.test.js
//
// A path to another build of the module can be given as the one argument.
//
// The buffers are written and read the way docs/assets/karyon-wasm.js writes
// and reads them for `render`: one length-prefixed buffer in, allocated with
// `alloc` and freed with `dealloc`, and one buffer out whose first byte says
// whether the rest is an answer or a message, freed by its own length.

const fs = require("fs");
const path = require("path");
const assert = require("assert");

const WASM =
  process.argv[2] ||
  path.join(__dirname, "..", "target", "wasm32-unknown-unknown", "release", "karyon_playground.wasm");

const encoder = new TextEncoder();
const decoder = new TextDecoder();

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

// ------------------------------------------------------------------ buffers

function pack(write) {
  const parts = [];
  let total = 0;
  const into = {
    u32(n) {
      const b = new Uint8Array(4);
      new DataView(b.buffer).setUint32(0, n, true);
      parts.push(b);
      total += 4;
    },
    str(s) {
      const bytes = encoder.encode(s);
      into.u32(bytes.length);
      parts.push(bytes);
      total += bytes.length;
    },
  };
  write(into);
  const out = new Uint8Array(total);
  let at = 0;
  parts.forEach((part) => {
    out.set(part, at);
    at += part.length;
  });
  return out;
}

// Calls one export with one buffer in, and hands back {ok, bytes} out.
function call(wasm, exported, input) {
  const inPtr = wasm.alloc(input.length);
  new Uint8Array(wasm.memory.buffer, inPtr, input.length).set(input);
  const outPtr = wasm[exported](inPtr, input.length);
  wasm.dealloc(inPtr, input.length);
  // Read again rather than reused from above: the buffer is replaced whenever
  // the module grows its memory, and drawing a figure is what grows it.
  const view = new DataView(wasm.memory.buffer);
  const ok = view.getUint8(outPtr) === 1;
  const len = view.getUint32(outPtr + 1, true);
  const bytes = new Uint8Array(wasm.memory.buffer.slice(outPtr + 5, outPtr + 5 + len));
  wasm.dealloc(outPtr, 5 + len);
  return { ok, bytes };
}

function figures(wasm) {
  const answer = call(wasm, "figures", new Uint8Array(0));
  assert.ok(answer.ok, "figures answered with a message");
  const view = new DataView(answer.bytes.buffer);
  let at = 0;
  const u32 = () => {
    const n = view.getUint32(at, true);
    at += 4;
    return n;
  };
  const str = () => {
    const len = u32();
    const text = decoder.decode(answer.bytes.subarray(at, at + len));
    at += len;
    return text;
  };
  const count = u32();
  const list = [];
  for (let i = 0; i < count; i++) {
    const name = str();
    const moves = answer.bytes[at] === 1;
    at += 1;
    const region = str();
    list.push({ name, moves, region });
  }
  assert.strictEqual(at, answer.bytes.length, "the list is exactly as long as it says");
  return list;
}

function figure(wasm, name, theme, background, width, region, prefix) {
  const input = pack((into) => {
    into.str(name);
    into.str(theme);
    into.str(background);
    into.u32(width);
    into.str(region);
    into.str(prefix);
  });
  const answer = call(wasm, "figure", input);
  return { ok: answer.ok, body: decoder.decode(answer.bytes) };
}

// ------------------------------------------------------------------- checks

if (!fs.existsSync(WASM)) {
  console.log("no module at " + WASM + "; build it first, as the comment at the top says");
  process.exit(1);
}

WebAssembly.instantiate(fs.readFileSync(WASM), {}).then((built) => {
  const wasm = built.instance.exports;
  console.log("module " + WASM + ", " + fs.statSync(WASM).size + " bytes");

  let started = process.hrtime.bigint();
  const list = figures(wasm);
  const listed = Number(process.hrtime.bigint() - started) / 1e6;
  console.log("figures: " + list.length + " in " + listed.toFixed(0) + " ms");
  for (const entry of list) {
    console.log("  " + entry.name.padEnd(36) + (entry.moves ? entry.region : "-"));
  }

  check("every committed figure is listed", () => {
    const committed = fs
      .readdirSync(path.join(__dirname, "..", "..", "assets"))
      .filter((file) => file.endsWith(".svg"))
      .map((file) => file.slice(0, -4))
      .sort();
    assert.deepStrictEqual(list.map((entry) => entry.name).sort(), committed);
  });

  started = process.hrtime.bigint();
  const drawn = figure(wasm, "example-genomewide", "dark", "", 700, "", "live-3-");
  const took = Number(process.hrtime.bigint() - started) / 1e6;
  console.log(
    "\nfigure example-genomewide, dark, 700 px, prefix live-3-, " + took.toFixed(0) + " ms:"
  );
  console.log(drawn.body.slice(0, 300));
  console.log("");

  check("it is a figure, 700 pixels wide, on the dark page", () => {
    assert.ok(drawn.ok, drawn.body);
    assert.ok(drawn.body.startsWith('<svg xmlns="http://www.w3.org/2000/svg" width="700" '));
    assert.ok(drawn.body.includes('fill="#14181d"'));
  });

  check("every id it writes carries the prefix, and it does not name itself", () => {
    const ids = [...drawn.body.matchAll(/ id="([^"]*)"/g)].map((m) => m[1]);
    assert.ok(ids.length > 0);
    assert.ok(ids.every((id) => id.startsWith("live-3-")), ids.join(" "));
    assert.ok(!drawn.body.slice(0, drawn.body.indexOf(">")).includes("aria-labelledby"));
  });

  check("asked for nothing, a figure is the committed file", () => {
    for (const name of ["example", "example-genomewide", "gallery", "example-maps"]) {
      const plain = figure(wasm, name, "light", "", 0, "", "");
      const file = fs.readFileSync(path.join(__dirname, "..", "..", "assets", name + ".svg"), "utf8");
      assert.ok(plain.ok, plain.body);
      assert.ok(plain.body === file, name + " differs from assets/" + name + ".svg");
    }
  });

  check("a window half as wide moves the ruler", () => {
    const whole = figure(wasm, "example", "light", "", 0, "", "");
    const half = figure(wasm, "example", "light", "", 0, "NC_000962.3:761000-761999", "");
    assert.ok(half.ok, half.body);
    assert.notStrictEqual(half.body, whole.body);
    assert.ok(half.body.includes(">NC_000962.3:761000-761999</text>"));
  });

  check("a request that makes no sense is a message, not a trap", () => {
    const wrong = figure(wasm, "example", "sepia", "", 0, "", "");
    assert.ok(!wrong.ok);
    assert.ok(wrong.body.includes("light or dark"), wrong.body);
  });

  if (failures) {
    console.log(failures + " failed");
    process.exit(1);
  }
});
