---
hide:
  - navigation
  - toc
---

<link rel="stylesheet" href="../stylesheets/tree-viewer.css">

# Tree viewer

<noscript>
  <p><strong>This page needs JavaScript.</strong> Everything it does happens in
  your browser, so there is nothing for it to fall back to. The same views come
  out of a terminal: see <a href="../guide/cli/">Command line</a>.</p>
</noscript>

<div class="tv-drop" id="tv-drop" hidden markdown="0">
  <p>Try it on a tree made here and now:
  <button class="tv-btn" type="button" data-example="200" data-label="200 tips">200 tips</button>
  <button class="tv-btn" type="button" data-example="20000" data-label="20,000 tips">20,000 tips</button>
  <button class="tv-btn" type="button" data-example="1000000" data-label="a million tips">a million tips</button>
  </p>
  <p>Drop a Newick file anywhere on this page, or
  <label class="tv-file">choose one<input type="file" id="tv-file" accept=".nwk,.newick,.tree,.tre,.txt"></label>.
  Nothing leaves your browser. Drop a second file with a header and one row per tip and it is drawn as
  strips beside the names.</p>
  <details>
    <summary>Or paste one</summary>
    <textarea id="tv-paste" rows="4" placeholder="((A:0.1,B:0.2):0.3,C:0.4);"></textarea>
    <button class="tv-btn" id="tv-usepaste" type="button">Draw it</button>
  </details>
</div>

<div class="tv-app" id="tv-app" hidden markdown="0">
  <div class="tv-bar">
    <button class="tv-btn" id="tv-shape" type="button" aria-pressed="false">Phylogram</button>
    <input class="tv-search" id="tv-search" type="search" placeholder="Find a tip, then Enter">
    <output class="tv-count" id="tv-count"></output>
    <output class="tv-rows" id="tv-rowsout"></output>
    <output class="tv-detail" id="tv-detail"></output>
    <span class="tv-sheet" id="tv-sheet" hidden>
      <span id="tv-sheetname"></span>
      <button class="tv-btn" id="tv-dropsheet" type="button" title="Draw without it">Drop</button>
    </span>
    <div class="tv-spacer"></div>
    <button class="tv-btn" id="tv-fit" type="button">Fit</button>
    <button class="tv-btn" id="tv-export" type="button" title="Save this view as karyon's own SVG">Export SVG</button>
  </div>

  <p class="tv-error" id="tv-error" hidden></p>

  <div class="tv-plot" id="tv-plot">
    <canvas id="tv-canvas"></canvas>
  </div>

  <p class="tv-hint">Drag to move, wheel to zoom, double-click to zoom in. The layout is worked out once by
    the program and never again: what moves is the window onto it, so a gesture costs a repaint of what is on
    screen and never a walk over the tree. Export hands the view back to the program and saves karyon's own
    figure of it.</p>

  <p class="tv-command">The same view from a shell: <code id="tv-command"></code></p>
</div>

<script src="../assets/karyon-wasm.js"></script>
<script src="../assets/tree-canvas.js"></script>
<script src="../assets/tree-viewer.js"></script>
