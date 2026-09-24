---
hide:
  - navigation
  - toc
---

<link rel="stylesheet" href="../stylesheets/playground.css">

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
      <div class="pg-plot" id="pg-plot" aria-live="polite"></div>
      <div class="pg-controls" id="pg-controls"></div>
      <div class="pg-foot"><span class="pg-status" id="pg-status" role="status" aria-live="polite">loading the program…</span></div>
    </section>
  </div>
</div>

<dialog class="pg-picker" id="pg-picker" aria-labelledby="pg-picker-title">
  <aside class="pg-panel">
    <div class="pg-panel-head">
      <div>
        <h2 id="pg-picker-title">Examples</h2>
        <p>Twenty-four examples that between them use all twenty-eight track flags. Every preview is drawn live, in this page.</p>
      </div>
      <button class="pg-chip" id="pg-picker-close" type="button">Close</button>
    </div>
    <input class="pg-search" id="pg-search" type="search" placeholder="Search examples" autocomplete="off" aria-label="Search examples">
    <div class="pg-panel-body" id="pg-picker-body"></div>
  </aside>
</dialog>

<div class="pg-notes" markdown>

## How to use it

1. **Start from an example.** *Examples* opens two dozen of them, which between
   them use every track the command line has. Or type your own command in the
   box: it takes everything you would type after `karyon` in a terminal.
2. **Edit the input files** in the tabs under the command. The command refers to
   a file by its tab name. `+` adds a file; double-click a tab (or press F2) to
   rename it.
3. **Press Draw.** Turn on *Interactive* to drag the figure along the genome and
   scroll to zoom once it has focus; every frame is a fresh run of the command.
4. **Use the controls under the figure.** They depend on the example, and each
   one rewrites its flag in the command, so the command always describes what
   you see.
5. **Export SVG** saves the figure exactly as drawn.

Everything runs in your browser: karyon is compiled to WebAssembly, and nothing
you type or paste is uploaded anywhere.

A few differences from a terminal:

- The page picks the figure width and follows the dark theme. Put `--width` or
  `--theme` in the command to override either.
- There is no standard input (`-`) and no `--help` here. The
  [Command line](guide/cli.md) guide is the help text written out.

</div>

<script src="../assets/karyon-wasm.js" defer></script>
<script src="../assets/playground.js" defer></script>
