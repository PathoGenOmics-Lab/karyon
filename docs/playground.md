---
hide:
  - navigation
  - toc
---

# Playground

<div id="karyon-playground" markdown="0">
  <noscript>
    <p><strong>This page needs JavaScript.</strong> Everything it does is done in
    your browser, so there is nothing for it to fall back to. The same commands
    run in a terminal: see <a href="../guide/cli/">Command line</a>.</p>
  </noscript>
</div>

<div class="pg-app" id="pg-app" hidden markdown="0">
  <div class="pg-bar">
    <div class="pg-group">
      <button class="pg-btn" id="pg-examples" type="button" aria-haspopup="dialog" aria-expanded="false">Examples</button>
    </div>
    <div class="pg-group">
      <button class="pg-btn pg-primary" id="pg-draw" type="button">Draw</button>
      <label class="pg-toggle"><input type="checkbox" id="pg-live"> Interactive</label>
    </div>
    <div class="pg-spacer"></div>
    <output class="pg-region" id="pg-region"></output>
    <div class="pg-group">
      <button class="pg-btn" id="pg-reset" type="button" title="Back to the region the command names">Reset view</button>
      <button class="pg-btn" id="pg-layout" type="button" title="Side by side or stacked">Layout</button>
      <button class="pg-btn" id="pg-export" type="button" title="Save the figure as it stands">Export SVG</button>
      <button class="pg-btn" id="pg-full" type="button" title="Fill the screen">Full screen</button>
    </div>
  </div>

  <div class="pg-panes" id="pg-panes">
    <section class="pg-editor" aria-label="Input">
      <label class="pg-legend" for="pg-command">Command</label>
      <textarea id="pg-command" class="pg-code pg-command" rows="4" spellcheck="false" aria-label="Command"></textarea>

      <div class="pg-tabs" id="pg-tabs" role="tablist" aria-label="Files"></div>
      <textarea id="pg-file" class="pg-code pg-body" spellcheck="false" aria-label="File contents"></textarea>
    </section>

    <div class="pg-split" id="pg-split" role="separator" aria-orientation="vertical" aria-valuenow="38" aria-valuemin="20" aria-valuemax="75" tabindex="0" aria-label="Resize the panes"></div>

    <section class="pg-view" aria-label="Figure">
      <div class="pg-plot" id="pg-plot" role="img" aria-label="The figure this command draws"></div>
      <div class="pg-foot"><span class="pg-status" id="pg-status" role="status" aria-live="polite">loading the program…</span></div>
    </section>
  </div>
</div>

<dialog class="pg-picker" id="pg-picker" aria-labelledby="pg-picker-title">
  <aside class="pg-panel">
    <div class="pg-panel-head">
      <div>
        <h2 id="pg-picker-title">Examples</h2>
        <p>Twenty-one, between them every flag the command has. Each preview is drawn by the program, here, as this opened.</p>
      </div>
      <button class="pg-chip" id="pg-picker-close" type="button">Close</button>
    </div>
    <input class="pg-search" id="pg-search" type="search" placeholder="Search examples" autocomplete="off" aria-label="Search examples">
    <div class="pg-panel-body" id="pg-picker-body"></div>
  </aside>
</dialog>

<div class="pg-notes" markdown>

Every flag the command has and every reader behind it, over the files on the
left. The region string, the coordinate conventions, the counts and every
refusal are what the terminal gives, because it is the same code: the grammar
lives in the library and takes a closure that answers with a file's text, so a
shell hands it a disk and this page hands it the editor.

**Interactive** re-runs the whole program on every frame. Drag the figure to
pan, scroll to zoom, and the command's region string follows along, because
that is what is actually being changed. Nothing is transformed or scaled: every
frame is a figure `karyon` drew, at that region, from those files.

The pane supplies two things the command does not: the width to draw at, and
the dark theme when the page is dark. Writing `--width` or `--theme` in the box
overrules both, since the command is the thing that decides.

Two things are a terminal's and are answered as text here. `-` for standard
input, because nothing is piped into a page. And `--help`, which is
[the command line guide](guide/cli.md) written out.

</div>

<script src="../assets/karyon-wasm.js" defer></script>
<script src="../assets/playground.js" defer></script>
