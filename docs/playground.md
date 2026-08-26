---
hide:
  - toc
---

# Playground

The command line, running in this page. Edit either box and press **Draw**.

Nothing is uploaded and nothing is downloaded past the first load: the whole of
`karyon` is compiled to WebAssembly and the figure is drawn by your own browser,
from the files in the box below. It works with the network unplugged.

<div id="karyon-playground" markdown="0">
  <noscript>
    <p><strong>This page needs JavaScript.</strong> Everything it does is done
    in your browser, so there is nothing for it to fall back to. The same
    commands run in a terminal: see <a href="../guide/cli/">Command line</a>.</p>
  </noscript>
</div>

<div id="karyon-playground-app" markdown="0" hidden>
  <div class="pg-examples">
    <span class="pg-label">Start from</span>
    <span id="pg-examples-buttons"></span>
  </div>

  <label class="pg-label" for="pg-command">Command</label>
  <textarea id="pg-command" class="pg-input" rows="3" spellcheck="false"></textarea>

  <label class="pg-label" for="pg-files">Files, one <code>=== name ===</code> header each</label>
  <textarea id="pg-files" class="pg-input" rows="10" spellcheck="false"></textarea>

  <div class="pg-actions">
    <button id="pg-draw" class="pg-button" type="button">Draw</button>
    <span id="pg-status" class="pg-status">loading the program…</span>
  </div>

  <div id="pg-output" class="pg-output"></div>
</div>

## What it can and cannot do

Every flag the command has, and every reader behind them, over the files in the
box. The region string, the coordinate conventions, the counts and every
refusal are the same code the terminal runs, because it is the same code: the
grammar lives in the library and takes a closure that answers with a file's
text, so a shell hands it a disk and this page hands it the box above.

Two things are a terminal's and are answered as text here instead. `-` for
standard input, because nothing is piped into a page. And `--help`, which is
[the command line guide](guide/cli.md) written out.

The figure below is a real `<svg>` element, so it is selectable, searchable and
readable by a screen reader, and **right click, save image as** keeps it. The
program is about half a megabyte and is fetched once.

<!-- Relative to this page, which MkDocs publishes at `playground/index.html`,
     so `../assets/` is the site's own assets directory whatever the site is
     served under. `assets/` alone would look inside the page's own folder. -->
<script src="../assets/playground.js" defer></script>
