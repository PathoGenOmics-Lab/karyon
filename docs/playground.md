---
hide:
  - navigation
  - toc
---

<link rel="stylesheet" href="../stylesheets/playground.css">

# Playground

Write a karyon command, edit the files it reads and watch the figure it draws,
with nothing to install. This page runs the same program as the command line,
compiled to WebAssembly.
{ .k-lead }

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
        <p>Twenty-four examples that between them use twenty-eight of the track flags. Every preview is drawn live, in this page.</p>
      </div>
      <button class="pg-chip" id="pg-picker-close" type="button">Close</button>
    </div>
    <input class="pg-search" id="pg-search" type="search" placeholder="Search examples" autocomplete="off" aria-label="Search examples">
    <div class="pg-panel-body" id="pg-picker-body"></div>
  </aside>
</dialog>

<div class="pg-notes" markdown>

## How to use it

1. **Start from an example, or write a command.** *Examples* opens twenty-four
   figures that between them use every track flag the command line has. The
   command box takes everything you would type after `karyon` in a terminal.
2. **Edit the files.** Each tab under the command is one input file, and the
   command names it by the tab's name. `+` adds a file, a double click or
   ++f2++ renames one, and ++delete++ or the cross on the tab removes it.
3. **Draw.** The figure redraws a moment after you stop typing, and *Draw* or
   ++ctrl+enter++ draws it at once.
4. **Explore.** Tick *Interactive*, then drag the figure along the sequence,
   or give it focus and scroll, press the arrow keys, or press ++plus++ and
   ++minus++. Every frame is a fresh run of the command, and the region in
   the command follows you. *Reset view* puts back the command the example
   started with.
5. **Use the controls under the figure.** They depend on the example. The
   window slider and the flag controls rewrite the command, so it always says
   what you see; the others regenerate the example's input files.
6. **Export SVG** saves the figure exactly as drawn. *Layout* puts the command
   beside the figure or above it, and *Full screen* gives it the whole screen.

!!! note "Nothing leaves your browser"
    The program runs inside this page, and nothing you type or paste is
    uploaded anywhere.

A few things differ from a terminal:

- **The page picks the width and the theme.** The figure is drawn as wide as
  its pane and follows the site's dark setting. Put `--width` or `--theme` in
  the command to override either.
- **Files are the tabs.** A file name in the command names a tab, not a file
  on your disk.
- **There is no standard input and no help text.** `-` has nothing to read
  here, and `--help` answers with a note. The
  [Command line](guide/cli.md) guide is the help text written out.

</div>

<script src="../assets/karyon-wasm.js" defer></script>
<script src="../assets/playground.js" defer></script>
