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
  Nothing leaves your browser.</p>
  <details>
    <summary>Or paste one</summary>
    <textarea id="tv-paste" rows="4" placeholder="((A:0.1,B:0.2):0.3,C:0.4);"></textarea>
    <button class="tv-btn" id="tv-usepaste" type="button">Draw it</button>
  </details>
</div>

<div class="tv-app" id="tv-app" hidden markdown="0">
  <div class="tv-bar">
    <div class="tv-group" role="group" aria-label="Projection">
      <button class="tv-btn" type="button" data-projection="rectangular" aria-pressed="true">Rectangular</button>
      <button class="tv-btn" type="button" data-projection="circular" aria-pressed="false">Circular</button>
      <button class="tv-btn" type="button" data-projection="unrooted" aria-pressed="false">Unrooted</button>
    </div>
    <label class="tv-rows">
      <input type="range" id="tv-rows" min="8" max="400" step="4" value="60">
      <output id="tv-rowsout">60 rows</output>
    </label>
    <input class="tv-search" id="tv-search" type="search" placeholder="Find a tip, then Enter">
    <output class="tv-timing" id="tv-timing"></output>
  </div>

  <nav class="tv-trail" id="tv-trail" aria-label="Where you are"></nav>

  <p class="tv-error" id="tv-error" hidden></p>

  <div class="tv-plot" id="tv-plot">
    <div class="tv-stage" id="tv-stage"></div>
    <div class="tv-hud">
      <output class="tv-zoom" id="tv-zoom"></output>
      <button class="tv-btn" id="tv-fit" type="button" disabled>Fit</button>
    </div>
  </div>
  <p class="tv-hint">Drag to move, wheel to zoom, double-click to zoom in, and keep pulling the wheel back
    to come out a level at a time. Stop moving and the tree is redrawn at the detail the view is asking
    for, so five gestures take twenty thousand tips down to five hundred with every name legible.</p>

  <p class="tv-command">The same view from a shell: <code id="tv-command"></code></p>
</div>

<script src="../assets/karyon-wasm.js"></script>
<script src="../assets/tree-viewer.js"></script>
