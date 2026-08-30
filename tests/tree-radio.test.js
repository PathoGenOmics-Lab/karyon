// The row of buttons where one of them is on, checked by running it. Plain
// node, no framework and no install, like the canvas beside it. Run with
// `node tests/tree-radio.test.js`.

const path = require("path");
const assert = require("assert");

// ---------------------------------------------------------- enough of a DOM

// Only what the control touches. Anything it reaches for that is not here is a
// dependency it should not have.
function fakeButton(attribute, value, checked) {
  const attributes = { [attribute]: value, "aria-checked": String(!!checked) };
  return {
    textContent: value,
    tabIndex: checked ? 0 : -1,
    focused: 0,
    getAttribute: (name) => (name in attributes ? attributes[name] : null),
    setAttribute(name, what) { attributes[name] = what; },
    focus() { this.focused += 1; },
    closest() { return this; },
  };
}

function fakeGroup(attribute, values, checkedAt) {
  const buttons = values.map((value, at) => fakeButton(attribute, value, at === checkedAt));
  const listeners = {};
  return {
    buttons,
    querySelectorAll: () => buttons,
    addEventListener(kind, run) { listeners[kind] = run; },
    click(button) { listeners.click({ target: button }); },
    press(key, on) {
      let stopped = 0;
      listeners.keydown({ key, target: on, preventDefault() { stopped += 1; } });
      return stopped;
    },
  };
}

global.window = {};
require(path.join(__dirname, "..", "docs", "assets", "tree-radio.js"));
const radio = global.window.karyonRadio;

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

const marks = (box) => box.buttons.map((b) => `${b.textContent}=${b.getAttribute("aria-checked")}/${b.tabIndex}`);

check("one of them is on, and it is the only one a tab reaches", () => {
  const box = fakeGroup("data-projection", ["rows", "disc", "spread"], 0);
  radio.make(box, "data-projection", () => {});
  assert.deepStrictEqual(marks(box), ["rows=true/0", "disc=false/-1", "spread=false/-1"]);
  box.click(box.buttons[2]);
  assert.deepStrictEqual(marks(box), ["rows=false/-1", "disc=false/-1", "spread=true/0"]);
});

check("and the check bites: leaving every button reachable fails it", () => {
  const box = fakeGroup("data-projection", ["rows", "disc", "spread"], 0);
  radio.make(box, "data-projection", () => {});
  box.click(box.buttons[1]);
  const reachable = box.buttons.filter((b) => b.tabIndex === 0);
  assert.strictEqual(reachable.length, 1, `${reachable.length} buttons are reachable by tab`);
});

check("taking the one already taken says nothing happened", () => {
  const box = fakeGroup("data-lengths", ["phylogram", "cladogram"], 0);
  const said = [];
  radio.make(box, "data-lengths", (value) => said.push(value));
  box.click(box.buttons[0]);
  assert.deepStrictEqual(said, [], "it reported a change that did not happen");
  box.click(box.buttons[1]);
  assert.deepStrictEqual(said, ["cladogram"]);
  box.click(box.buttons[1]);
  assert.deepStrictEqual(said, ["cladogram"], "it reported the same change twice");
});

check("the arrow keys move through them and take what they land on", () => {
  const box = fakeGroup("data-projection", ["rows", "disc", "spread"], 0);
  const said = [];
  radio.make(box, "data-projection", (value) => said.push(value));
  box.press("ArrowRight", box.buttons[0]);
  box.press("ArrowRight", box.buttons[1]);
  assert.deepStrictEqual(said, ["disc", "spread"]);
  assert.strictEqual(box.buttons[2].focused, 1, "the key moved the mark but not the focus");
  // Round the end and back to the start.
  box.press("ArrowRight", box.buttons[2]);
  assert.deepStrictEqual(said, ["disc", "spread", "rows"], "it did not come round");
  box.press("ArrowLeft", box.buttons[0]);
  assert.deepStrictEqual(said, ["disc", "spread", "rows", "spread"], "it did not go back round");
  box.press("Home", box.buttons[2]);
  box.press("End", box.buttons[0]);
  assert.deepStrictEqual(marks(box), ["rows=false/-1", "disc=false/-1", "spread=true/0"]);
});

check("a key it has no use for is left for the page", () => {
  const box = fakeGroup("data-projection", ["rows", "disc"], 0);
  const said = [];
  radio.make(box, "data-projection", (value) => said.push(value));
  assert.strictEqual(box.press("Tab", box.buttons[0]), 0, "it swallowed Tab");
  assert.strictEqual(box.press("a", box.buttons[0]), 0, "it swallowed a letter");
  assert.strictEqual(box.press("ArrowRight", box.buttons[0]), 1, "it did not take the key it uses");
  assert.deepStrictEqual(said, ["disc"]);
});

check("a key pressed somewhere else in the group is not its business", () => {
  const box = fakeGroup("data-projection", ["rows", "disc"], 0);
  const said = [];
  radio.make(box, "data-projection", (value) => said.push(value));
  assert.strictEqual(box.press("ArrowRight", { textContent: "not a button" }), 0);
  assert.deepStrictEqual(said, [], "it acted on a key pressed on something else");
});

process.exit(failures ? 1 : 0);
