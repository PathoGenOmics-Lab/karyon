// A row of buttons where one of them is on.
//
// The browser has no such control. A single button that cycles through the
// states cannot say whether its label is where you are or where you would go,
// and this page had two of them. What it has instead is what a reader who knows
// radio buttons will already try: the one that is on is the only one a tab
// reaches, and the arrow keys move between them and take the one they land on.
//
// It is here rather than inside the viewer so that it can be run without a
// browser, which is the only way anything on this page gets checked.

window.karyonRadio = (function () {
  "use strict";

  // `box` holds the buttons, each carrying `attribute` with the value it
  // stands for. `pick` is called with that value whenever one is taken, and not
  // when the one already taken is taken again.
  function make(box, attribute, pick) {
    var buttons = [].slice.call(box.querySelectorAll("[" + attribute + "]"));
    var taken = null;

    function show(value) {
      taken = value;
      for (var i = 0; i < buttons.length; i++) {
        var on = buttons[i].getAttribute(attribute) === value;
        buttons[i].setAttribute("aria-checked", String(on));
        buttons[i].tabIndex = on ? 0 : -1;
      }
    }

    function take(button, focus) {
      if (!button) return;
      var value = button.getAttribute(attribute);
      var changed = value !== taken;
      show(value);
      if (focus && button.focus) button.focus();
      if (changed) pick(value);
    }

    for (var at = 0; at < buttons.length; at++) {
      if (buttons[at].getAttribute("aria-checked") === "true") {
        taken = buttons[at].getAttribute(attribute);
      }
    }

    box.addEventListener("click", function (event) {
      var button = event.target && event.target.closest
        ? event.target.closest("[" + attribute + "]")
        : null;
      take(button, false);
    });

    box.addEventListener("keydown", function (event) {
      var at = buttons.indexOf(event.target);
      if (at < 0) return;
      var to = at;
      if (event.key === "ArrowRight" || event.key === "ArrowDown") to = at + 1;
      else if (event.key === "ArrowLeft" || event.key === "ArrowUp") to = at - 1;
      else if (event.key === "Home") to = 0;
      else if (event.key === "End") to = buttons.length - 1;
      else return;
      if (event.preventDefault) event.preventDefault();
      take(buttons[(to + buttons.length) % buttons.length], true);
    });

    return {
      show: show,
      taken: function () { return taken; },
    };
  }

  return { make: make };
})();
