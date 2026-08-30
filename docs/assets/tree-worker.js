// The program, off the page's thread.
//
// Drawing a tree is not free: a million tip figure takes a couple of seconds
// even with the file already read, and on the page's own thread those seconds
// are seconds in which nothing moves, no wheel turns and no drag lands. So the
// program runs here instead. The page keeps the hand, which is a transform and
// costs a composite; the worker keeps the drawing.
//
// The files are sent once and kept. A tree is megabytes of text and posting it
// with every request would put the copying back on the page's thread, which is
// the cost this exists to remove; keeping it here also means the program's own
// memory of the last tree it read hits on every call rather than never.

/* global importScripts, karyon */
importScripts("karyon-wasm.js");

var files = [];

function answer(message) {
  self.postMessage(message);
}

self.onmessage = function (event) {
  var job = event.data;

  if (job.kind === "files") {
    files = job.files || [];
    return;
  }

  karyon.load().then(function () {
    if (job.kind === "draw") {
      var drawn = karyon.run(job.command, files, job.room);
      answer({ kind: "drawn", id: job.id, ok: drawn.ok, body: drawn.body, ms: drawn.ms });
      return;
    }

    if (job.kind === "choose") {
      // The candidates a zoom is deciding between, widest first, answered with
      // the first that draws less of the tree than is already on screen. Done
      // here and not on the page because it is up to three drawings, and three
      // drawings of a million tips is six seconds of a page that cannot be
      // touched.
      for (var i = 0; i < job.candidates.length; i++) {
        var candidate = job.candidates[i];
        var drawn = karyon.run(candidate.command, files, job.room);
        if (!drawn.ok) continue;
        if (karyon.tipsAccountedFor(drawn.body) >= job.was) continue;
        answer({
          kind: "chosen",
          id: job.id,
          focus: candidate.focus,
          label: candidate.label,
          body: drawn.body,
          ms: drawn.ms,
        });
        return;
      }
      // Every one of them draws what is already there, so the view has nothing
      // below it to open.
      answer({ kind: "chosen", id: job.id, focus: null });
    }
  });
};
