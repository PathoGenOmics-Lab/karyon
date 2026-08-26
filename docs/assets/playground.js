// The page half of the playground's bridge to WebAssembly.
//
// The protocol is written down in playground/src/lib.rs and repeated here only
// where the code would otherwise be a row of magic offsets: one buffer in and
// one buffer out, every number a little-endian u32, every string UTF-8.
//
// No framework and no build step, because the site is Material for MkDocs and
// the crate has no dependencies, and a playground that needed a bundler to
// demonstrate a program that needs nothing would be making the wrong point.

(function () {
  "use strict";

  // The protocol, the command line parsing and the region arithmetic are in
  // karyon-wasm.js, which the home page runs its own figure over too. Two
  // copies of a protocol are two things that drift.
  var K = window.karyon;

  var el = {};
  var files = [];
  var active = 0;
  var pending = null;

  var EXAMPLES = [
    {
      name: "A locus",
      command:
        "NC_000962.3:761,000-762,999 --coverage depth.bg --label depth --aggregate min \\\n" +
        "  --features genes.gff3 --label annotation \\\n" +
        "  --variants calls.vcf --label variants \\\n" +
        "  --title 'rpoB locus, resistance determining region'",
      // Wide enough that the window can be moved. The three rows that fitted
      // the region exactly refused on the first drag, which is the program
      // being right and the example being too small to show anything else.
      files: [
        {
          name: "depth.bg",
          body:
            "NC_000962.3 756999 759999 62\n" +
            "NC_000962.3 759999 760999 58\n" +
            "NC_000962.3 760999 761899 57\n" +
            "NC_000962.3 761899 762029 3\n" +
            "NC_000962.3 762029 763999 60\n" +
            "NC_000962.3 763999 766999 54\n",
        },
        {
          name: "genes.gff3",
          body:
            "##gff-version 3\n" +
            "NC_000962.3 . gene 759807 763325 . + . Name=rpoB\n" +
            "NC_000962.3 . gene 761082 761162 . + . Name=RRDR\n",
        },
        {
          name: "calls.vcf",
          body:
            "NC_000962.3 760106 . C T . . AF=0.09;ANN=T|synonymous_variant|LOW|rpoB\n" +
            "NC_000962.3 761052 . C T . . AF=0.12;ANN=T|synonymous_variant|LOW|rpoB\n" +
            "NC_000962.3 761109 . G T . . AF=0.98;ANN=T|missense_variant|MODERATE|rpoB\n" +
            "NC_000962.3 761139 . C T . . AF=0.55;ANN=T|missense_variant|MODERATE|rpoB\n" +
            "NC_000962.3 761155 . T C . . AF=1.00;ANN=C|missense_variant|MODERATE|rpoB\n" +
            "NC_000962.3 761156 . C T . . AF=0.21;ANN=T|synonymous_variant|LOW|rpoB\n" +
            "NC_000962.3 761606 . G A . . AF=0.07;ANN=A|synonymous_variant|LOW|rpoB\n" +
            "NC_000962.3 762206 . C T . . AF=0.15;ANN=T|synonymous_variant|LOW|rpoB\n",
        },
      ],
    },
    {
      name: "A whole chromosome",
      command:
        "chr1:1-2,000,000 --coverage depth.bg --label depth --aggregate min \\\n" +
        "  --windows gc.bg --label 'GC content' --style steps",
      files: [
        { name: "depth.bg", body: "# generated below\n" },
        { name: "gc.bg", body: "# generated below\n" },
      ],
      // Two thousand rows written here rather than typed out, so the example
      // is one a reader can pan and zoom rather than one they can read.
      make: function () {
        var depth = "";
        var gc = "";
        for (var i = 0; i < 2000; i++) {
          var at = i * 1000;
          var dip = i > 900 && i < 1000 ? 0.25 : 1;
          depth +=
            "chr1\t" + at + "\t" + (at + 1000) + "\t" +
            Math.round((40 + 18 * Math.sin(i / 40)) * dip) + "\n";
          gc +=
            "chr1\t" + at + "\t" + (at + 1000) + "\t" +
            (0.42 + 0.11 * Math.sin(i / 17)).toFixed(3) + "\n";
        }
        return [
          { name: "depth.bg", body: depth },
          { name: "gc.bg", body: gc },
        ];
      },
    },
    {
      name: "Two trees",
      command: "tangle:1-1000 --no-axis --tanglegram before.nwk --against after.nwk",
      files: [
        { name: "before.nwk", body: "((a:1,b:1):1,(c:1,d:1):1);\n" },
        { name: "after.nwk", body: "((a:1,c:1):1,(b:1,d:1):1);\n" },
      ],
    },
    {
      name: "Gene neighbourhoods",
      command: "ESX-1:1-4,000 --loci loci.bed --links hits.tsv --label 'ESX-1'",
      files: [
        {
          name: "loci.bed",
          body:
            "H37Rv\t0\t1200\tespA\t0\t+\n" +
            "H37Rv\t1300\t2100\tespC\t0\t+\n" +
            "H37Rv\t2200\t3000\tespD\t0\t-\n" +
            "CDC1551\t0\t1200\tespA2\t0\t+\n" +
            "CDC1551\t2200\t3000\tespD2\t0\t-\n" +
            "Erdman\t0\t1200\tespA3\t0\t+\n" +
            "Erdman\t1300\t2100\tespC3\t0\t+\n",
        },
        {
          name: "hits.tsv",
          body: "espA\tespA2\t99.1\nespD\tespD2\t97.4\nespA2\tespA3\t96.2\n",
        },
      ],
    },
    {
      name: "Protein domains",
      command: "protein:1-700 --domains domains.tsv --analysis Pfam",
      files: [
        {
          name: "domains.tsv",
          body:
            "PknB\tmd5\t626\tPfam\tPF00069\tProtein kinase domain\t11\t275\t1e-40\tT\t01-01-2026\n" +
            "PknB\tmd5\t626\tPfam\tPF03793\tPASTA domain\t341\t400\t3e-10\tT\t01-01-2026\n" +
            "PknB\tmd5\t626\tPfam\tPF03793\tPASTA domain\t410\t468\t4e-10\tT\t01-01-2026\n" +
            "PknD\tmd5\t664\tPfam\tPF00069\tProtein kinase domain\t14\t277\t9e-40\tT\t01-01-2026\n" +
            "PknE\tmd5\t565\tPfam\tPF00069\tProtein kinase domain\t16\t280\t3e-39\tT\t01-01-2026\n" +
            "PknE\tmd5\t565\tPfam\tPF03793\tPASTA domain\t400\t458\t7e-09\tT\t01-01-2026\n",
        },
      ],
    },
    {
      name: "One molecule at a time",
      command: "chr11:1-200 --bisulfite calls.txt --context CpG --label 'H19 ICR'",
      files: [
        {
          name: "calls.txt",
          body: (function () {
            var out = "";
            var sites = [12, 31, 48, 66, 89, 104, 127, 151, 168];
            for (var r = 1; r <= 10; r++) {
              var on = r <= 5;
              for (var s = 0; s < sites.length; s++) {
                if (r === 6 && s < 3) continue;
                var call = on ? "Z" : "z";
                out += "read" + r + "/1\t" + (on ? "+" : "-") + "\tchr11\t" +
                  (sites[s] + 1) + "\t" + call + "\n";
              }
            }
            return out;
          })(),
        },
      ],
    },
  ];

  // ---------------------------------------------------------------------
  // Drawing
  // ---------------------------------------------------------------------

  var drawn = null;

  function draw() {
    if (!K.ready()) return;
    save();
    var answer = K.run(el.command.value, files, el.plot.clientWidth - 24);

    if (answer.ok) {
      el.plot.className = "pg-plot" + (el.live.checked ? " pg-live" : "");
      el.plot.innerHTML = answer.body;
      drawn = answer.body;
      var where = K.locus(el.command.value);
      el.region.textContent = where
        ? where.seq + ":" + K.grouped(where.start) + "-" + K.grouped(where.end) +
          "  (" + K.grouped(where.end - where.start + 1) + " bases)"
        : "";
      el.status.textContent =
        files.length + (files.length === 1 ? " file" : " files") +
        ", drawn in " + answer.ms.toFixed(answer.ms < 10 ? 1 : 0) + " ms";
      el.status.classList.remove("pg-bad");
    } else {
      el.plot.className = "pg-plot pg-failed";
      el.plot.textContent = "karyon: " + answer.body;
      drawn = null;
      el.region.textContent = "";
      el.status.textContent = "the command did not draw";
      el.status.classList.add("pg-bad");
    }
  }

  function soon() {
    clearTimeout(pending);
    pending = setTimeout(draw, 200);
  }

  // ---------------------------------------------------------------------
  // Interaction: every frame is a figure the program drew
  // ---------------------------------------------------------------------

  var origin = null;

  function move(next) {
    el.command.value = next;
    draw();
  }

  function interactive(on) {
    el.plot.classList.toggle("pg-live", on);
    el.reset.disabled = !on;
  }

  function onDown(event) {
    if (!el.live.checked || event.button !== 0) return;
    origin = { x: event.clientX };
    el.plot.setPointerCapture(event.pointerId);
    el.plot.classList.add("pg-dragging");
  }

  function onMove(event) {
    if (!origin) return;
    var moved = event.clientX - origin.x;
    if (Math.abs(moved) < 1) return;
    origin.x = event.clientX;
    move(K.panned(el.command.value, moved / Math.max(1, el.plot.clientWidth)));
  }

  function onUp(event) {
    if (!origin) return;
    origin = null;
    el.plot.classList.remove("pg-dragging");
    if (el.plot.hasPointerCapture(event.pointerId)) {
      el.plot.releasePointerCapture(event.pointerId);
    }
  }

  function onWheel(event) {
    if (!el.live.checked) return;
    event.preventDefault();
    var box = el.plot.getBoundingClientRect();
    var at = Math.min(1, Math.max(0, (event.clientX - box.left) / box.width));
    move(K.zoomed(el.command.value, event.deltaY > 0 ? 1.25 : 0.8, at));
  }

  // ---------------------------------------------------------------------
  // Files, as tabs
  // ---------------------------------------------------------------------

  // Which file the editor is holding, which is not always the selected one:
  // it is nothing at all until a file has been put in it. Saving without that
  // distinction writes an empty editor over the first file of every example
  // the moment it is loaded.
  var showing = -1;

  function save() {
    if (showing >= 0 && files[showing]) files[showing].body = el.file.value;
  }

  function show(index) {
    save();
    active = Math.max(0, Math.min(index, files.length - 1));
    el.file.value = files.length ? files[active].body : "";
    el.file.disabled = !files.length;
    showing = files.length ? active : -1;
    tabs();
  }

  function tabs() {
    el.tabs.textContent = "";
    files.forEach(function (file, index) {
      var tab = document.createElement("button");
      tab.type = "button";
      tab.className = "pg-tab" + (index === active ? " pg-on" : "");
      tab.setAttribute("role", "tab");
      tab.setAttribute("aria-selected", index === active ? "true" : "false");
      tab.textContent = file.name;
      tab.addEventListener("click", function () { show(index); });
      // A double click renames it, since the name is what the command calls it
      // by and a file nobody can rename is a file with the wrong name.
      tab.addEventListener("dblclick", function () {
        var name = prompt("Name of this file, as the command calls it", file.name);
        if (name) { file.name = name.trim(); tabs(); draw(); }
      });
      var shut = document.createElement("span");
      shut.className = "pg-shut";
      shut.textContent = "×";
      shut.title = "Remove " + file.name;
      shut.addEventListener("click", function (event) {
        event.stopPropagation();
        files.splice(index, 1);
        show(Math.min(active, files.length - 1));
        draw();
      });
      tab.appendChild(shut);
      el.tabs.appendChild(tab);
    });

    var add = document.createElement("button");
    add.type = "button";
    add.className = "pg-tab pg-add";
    add.textContent = "+";
    add.title = "Add a file";
    add.addEventListener("click", function () {
      var name = prompt("Name of the new file, as the command will call it", "data.bed");
      if (!name) return;
      files.push({ name: name.trim(), body: "" });
      show(files.length - 1);
      el.file.focus();
    });
    el.tabs.appendChild(add);
  }

  var home = null;

  function load(example) {
    el.command.value = example.command;
    // The editor is holding the last example's file, not this one's, so it
    // has nothing to save.
    showing = -1;
    files = (example.make ? example.make() : example.files).map(function (file) {
      return { name: file.name, body: file.body };
    });
    active = 0;
    show(0);
    // `Reset view` means the region this example named, not the first one's.
    home = example.command;
    draw();
  }

  // ---------------------------------------------------------------------
  // Chrome
  // ---------------------------------------------------------------------

  function menu() {
    el.menu.textContent = "";
    EXAMPLES.forEach(function (example) {
      var item = document.createElement("button");
      item.type = "button";
      item.className = "pg-item";
      item.textContent = example.name;
      item.addEventListener("click", function () {
        el.menu.hidden = true;
        el.examples.setAttribute("aria-expanded", "false");
        load(example);
      });
      el.menu.appendChild(item);
    });
  }

  function exportSvg() {
    if (!drawn) return;
    var blob = new Blob([drawn], { type: "image/svg+xml" });
    var url = URL.createObjectURL(blob);
    var link = document.createElement("a");
    link.href = url;
    link.download = "karyon.svg";
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  function dragSplit() {
    var moving = false;
    el.split.addEventListener("pointerdown", function (event) {
      moving = true;
      el.split.setPointerCapture(event.pointerId);
    });
    el.split.addEventListener("pointermove", function (event) {
      if (!moving) return;
      var box = el.panes.getBoundingClientRect();
      var stacked = el.panes.classList.contains("pg-stacked");
      var fraction = stacked
        ? (event.clientY - box.top) / box.height
        : (event.clientX - box.left) / box.width;
      fraction = Math.min(0.75, Math.max(0.2, fraction));
      el.panes.style.setProperty("--pg-split", (fraction * 100).toFixed(1) + "%");
      draw();
    });
    el.split.addEventListener("pointerup", function (event) {
      moving = false;
      if (el.split.hasPointerCapture(event.pointerId)) {
        el.split.releasePointerCapture(event.pointerId);
      }
    });
  }

  function start() {
    var ids = ["app", "bar", "panes", "split", "command", "file", "tabs", "plot",
               "status", "region", "draw", "live", "reset", "layout", "export",
               "full", "examples", "menu"];
    ids.forEach(function (name) { el[name] = document.getElementById("pg-" + name); });
    if (!el.app || !el.command) return;
    el.app.hidden = false;

    menu();
    el.examples.addEventListener("click", function () {
      var open = el.menu.hidden;
      el.menu.hidden = !open;
      el.examples.setAttribute("aria-expanded", String(open));
    });
    document.addEventListener("click", function (event) {
      if (!el.menu.hidden && !el.menu.contains(event.target) && event.target !== el.examples) {
        el.menu.hidden = true;
        el.examples.setAttribute("aria-expanded", "false");
      }
    });

    el.draw.addEventListener("click", draw);
    el.command.addEventListener("input", soon);
    el.file.addEventListener("input", soon);
    [el.command, el.file].forEach(function (box) {
      box.addEventListener("keydown", function (event) {
        if ((event.ctrlKey || event.metaKey) && event.key === "Enter") {
          event.preventDefault();
          draw();
        }
      });
    });

    el.live.addEventListener("change", function () { interactive(el.live.checked); });
    el.plot.addEventListener("pointerdown", onDown);
    el.plot.addEventListener("pointermove", onMove);
    el.plot.addEventListener("pointerup", onUp);
    el.plot.addEventListener("pointercancel", onUp);
    el.plot.addEventListener("wheel", onWheel, { passive: false });

    el.reset.addEventListener("click", function () {
      if (home) { el.command.value = home; draw(); }
    });
    el.layout.addEventListener("click", function () {
      el.panes.classList.toggle("pg-stacked");
      draw();
    });
    el.export.addEventListener("click", exportSvg);
    el.full.addEventListener("click", function () {
      if (document.fullscreenElement) document.exitFullscreen();
      else if (el.app.requestFullscreen) el.app.requestFullscreen();
    });
    dragSplit();

    var wide = null;
    window.addEventListener("resize", function () {
      clearTimeout(wide);
      wide = setTimeout(draw, 150);
    });
    // The palette toggle rewrites an attribute rather than reloading, so the
    // figure has to be told.
    K.onScheme(draw);

    el.draw.disabled = true;
    interactive(false);
    load(EXAMPLES[0]);

    K.load()
      .then(function () {
        el.draw.disabled = false;
        draw();
      })
      .catch(function (error) {
        el.status.textContent = "the program did not load: " + error.message;
        el.status.classList.add("pg-bad");
        el.plot.className = "pg-plot pg-failed";
        el.plot.textContent =
          "This page needs assets/karyon_playground.wasm, which is built and " +
          "published with the site. The same commands run in a terminal.";
      });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
