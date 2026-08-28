/* Every figure on this site is drawn wider than the column it lands in. The
   widest is a gallery sheet of 3472 pixels shown at 624, which is 18 per cent,
   and at that size none of the labels inside it can be read. The figures are
   what the project has to show, so each one has to be openable at the size it
   was drawn.

   The control is the figure itself. It becomes a button, so it is reachable by
   keyboard and announced as something that can be activated, and what it opens
   is a native `<dialog>`, so the focus trap and the Escape key are the
   browser's rather than ours. Inside the dialog the figure starts fitted to the
   viewport and can be switched to its own size, which is the point: fitted is
   still larger than the column, and actual size is where the labels resolve.

   No dependencies, in a crate that has none. */

(function () {
  "use strict";

  var FIT = "fit";
  var FULL = "full";

  /* A figure worth opening is a committed karyon, and every one of those lives
     under assets/figures. The logo and the rest of the site's furniture do not,
     which is the test, rather than a list of the classes they happen to carry.
     Card thumbnails are skipped because they are already links to somewhere
     better, and the live figure because the program redraws it. */
  function candidates() {
    var found = [];
    var imgs = document.querySelectorAll(".md-typeset img");
    for (var i = 0; i < imgs.length; i++) {
      var img = imgs[i];
      var src = img.getAttribute("src");
      if (!src || src.indexOf("assets/figures/") === -1) continue;
      if (img.closest("a") || img.closest("button")) continue;
      if (img.closest(".k-live")) continue;
      found.push(img);
    }
    return found;
  }

  function build() {
    var dialog = document.createElement("dialog");
    dialog.className = "k-viewer";
    dialog.innerHTML =
      '<div class="k-viewer__bar">' +
      '<span class="k-viewer__name"></span>' +
      '<button type="button" class="k-viewer__size" aria-pressed="false"></button>' +
      '<a class="k-viewer__open" target="_blank" rel="noopener">Open the SVG</a>' +
      '<button type="button" class="k-viewer__shut" aria-label="Close">&times;</button>' +
      "</div>" +
      '<div class="k-viewer__stage" tabindex="0" role="group" ' +
      'aria-label="The figure. Use the arrow keys to move around it at actual size.">' +
      '<img alt=""></div>';
    document.body.appendChild(dialog);
    return dialog;
  }

  function ready() {
    var figures = candidates();
    if (!figures.length) return;

    var dialog = build();
    var stage = dialog.querySelector(".k-viewer__stage");
    var shown = dialog.querySelector(".k-viewer__stage img");
    var name = dialog.querySelector(".k-viewer__name");
    var size = dialog.querySelector(".k-viewer__size");
    var open = dialog.querySelector(".k-viewer__open");
    var mode = FIT;
    var opener = null;

    function setMode(next) {
      mode = next;
      stage.dataset.mode = mode;
      size.textContent = mode === FIT ? "Actual size" : "Fit to window";
      size.setAttribute("aria-pressed", mode === FULL ? "true" : "false");
      /* At actual size the picture is larger than the window and the stage is
         what scrolls, so that is where the arrow keys have to land. Chrome does
         not make a scrolling box focusable on its own. */
      if (mode === FULL) stage.focus();
    }

    function show(img) {
      opener = img;
      shown.src = img.currentSrc || img.src;
      shown.alt = img.alt || "";
      name.textContent = img.alt || "";
      open.href = img.currentSrc || img.src;
      setMode(FIT);
      stage.scrollTop = 0;
      stage.scrollLeft = 0;
      dialog.showModal();
    }

    for (var i = 0; i < figures.length; i++) {
      (function (img) {
        var button = document.createElement("button");
        button.type = "button";
        button.className = "k-zoom";
        button.setAttribute(
          "aria-label",
          img.alt ? "See at full size: " + img.alt : "See this figure at full size"
        );
        /* One committed figure draws its own dark page, so it must not be
           given the white card the others need. The project names it for what
           it demonstrates, and that name is the test. */
        if (/-dark\.svg$/.test(img.getAttribute("src"))) {
          button.dataset.ground = "dark";
        }
        img.parentNode.insertBefore(button, img);
        button.appendChild(img);
        button.addEventListener("click", function () {
          show(img);
        });
      })(figures[i]);
    }

    size.addEventListener("click", function () {
      setMode(mode === FIT ? FULL : FIT);
    });
    dialog.querySelector(".k-viewer__shut").addEventListener("click", function () {
      dialog.close();
    });

    /* Clicking the ground closes it. The stage fills the dialog, so the test is
       whether the pointer landed on the picture rather than beside it. */
    dialog.addEventListener("click", function (event) {
      if (event.target === shown || event.target.closest(".k-viewer__bar")) return;
      dialog.close();
    });

    dialog.addEventListener("close", function () {
      shown.removeAttribute("src");
      if (opener) opener.parentNode.focus();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", ready);
  } else {
    ready();
  }
})();
