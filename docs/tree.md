---
hide:
  - navigation
  - toc
---

<link rel="stylesheet" href="../stylesheets/tree-viewer.css">

# Tree viewer

Open a phylogenetic tree and explore it in your browser: rectangular, circular
or unrooted, from a handful of tips to a million. karyon lays the tree out once
and then only redraws what is on screen, so even very large trees stay smooth.
Nothing you open is uploaded.
{ .tv-intro }

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
    <textarea id="tv-paste" rows="4" aria-label="A tree in Newick"
      placeholder="((A:0.1,B:0.2):0.3,C:0.4);"></textarea>
    <button class="tv-btn" id="tv-usepaste" type="button">Draw it</button>
  </details>
</div>

<p class="tv-error" id="tv-error" role="alert" hidden markdown="0"></p>

<div class="tv-app" id="tv-app" hidden markdown="0">
  <div class="tv-bar">
    <div class="tv-group" id="tv-lengths" role="radiogroup" aria-label="Branch lengths">
      <button class="tv-btn" type="button" role="radio" aria-checked="true" tabindex="0" data-lengths="phylogram"
        title="Branches are drawn at the length the file gives them">Phylogram</button>
      <button class="tv-btn" type="button" role="radio" aria-checked="false" tabindex="-1" data-lengths="cladogram"
        title="Branches are counted rather than measured, and every tip lines up">Cladogram</button>
    </div>
    <div class="tv-group" id="tv-projection" role="radiogroup" aria-label="Projection">
      <button class="tv-btn" type="button" role="radio" aria-checked="true" tabindex="0" data-projection="rows"
        title="Root on the left, tips in rows on the right">Rectangular</button>
      <button class="tv-btn" type="button" role="radio" aria-checked="false" tabindex="-1" data-projection="disc"
        title="The same rows bent round a circle, with depth becoming radius">Circular</button>
      <button class="tv-btn" type="button" role="radio" aria-checked="false" tabindex="-1" data-projection="spread"
        title="No root and no rows: the topology laid out in a plane">Unrooted</button>
    </div>
    <input class="tv-search" id="tv-search" type="search" aria-label="Find tips by name; a part of a name works, and so does a list"
      placeholder="Find tips, then Enter">
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
    <label class="tv-btn tv-another" title="Read another Newick file">Another tree<input type="file" id="tv-another"
      accept=".nwk,.newick,.tree,.tre,.txt"></label>
  </div>

  <div class="tv-plot" id="tv-plot">
    <output class="tv-tip" id="tv-tip" hidden></output>
    <canvas id="tv-canvas" tabindex="0" role="img"
      aria-label="The phylogeny. Up and down move it, plus and minus zoom, Home fits it. In the circular and unrooted views the other two arrows move it as well."></canvas>
  </div>

  <ul class="tv-hint">
    <li><strong>Move:</strong> drag the tree, or tab to it and use the arrow keys.</li>
    <li><strong>Zoom:</strong> scroll, double-click, or press <kbd>+</kbd> and <kbd>-</kbd>. <strong>Fit</strong> or <kbd>Home</kbd> shows the whole tree again.</li>
    <li><strong>Find:</strong> type part of a name, or a list of names, in the search box and press <kbd>Enter</kbd>.</li>
    <li class="tv-map-note"><strong>Overview:</strong> on a wide screen a small picture of the whole tree marks the part you are looking at; click or drag it to jump there.</li>
    <li><strong>Traits:</strong> drop a table with a header and one row per tip to draw it as coloured strips beside the names.</li>
    <li><strong>Export SVG</strong> saves the current view as a karyon figure, and the line below gives the same view as a command.</li>
  </ul>

  <div class="tv-key" id="tv-key" hidden></div>

  <p class="tv-note" id="tv-note" role="status" hidden></p>

  <p class="tv-command">The same view from a shell: <code id="tv-command"></code></p>
</div>

<script src="../assets/karyon-wasm.js"></script>
<script src="../assets/tree-canvas.js"></script>
<script src="../assets/tree-radio.js"></script>
<script src="../assets/tree-viewer.js"></script>
