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

  // Worked out from where this script itself was loaded from, rather than
  // written down. The site is served from a sub-path on GitHub Pages and from
  // the root on a laptop, and a page under a folder resolves a bare relative
  // path against the folder rather than against the site.
  var here = (document.currentScript && document.currentScript.src) || "";
  var WASM = here
    ? here.replace(/[^/]*$/, "karyon_playground.wasm")
    : "../assets/karyon_playground.wasm";
  var wasm = null;
  var encoder = new TextEncoder();
  var decoder = new TextDecoder();

  var EXAMPLES = [
    {
      name: "A locus",
      command:
        "NC_000962.3:761,000-761,200 --coverage depth.bg --label depth \\\n" +
        "  --features genes.gff3 --label annotation \\\n" +
        "  --variants calls.vcf --label variants \\\n" +
        "  --title 'rpoB resistance determining region'",
      files:
        "=== depth.bg ===\n" +
        "NC_000962.3\t760999\t761080\t54\n" +
        "NC_000962.3\t761080\t761120\t11\n" +
        "NC_000962.3\t761120\t761200\t58\n" +
        "\n" +
        "=== genes.gff3 ===\n" +
        "##gff-version 3\n" +
        "NC_000962.3\t.\tgene\t761040\t761160\t.\t+\t.\tName=RRDR\n" +
        "\n" +
        "=== calls.vcf ===\n" +
        "NC_000962.3\t761110\t.\tC\tT\t.\t.\t.\n" +
        "NC_000962.3\t761155\t.\tG\tA\t.\t.\t.\n",
    },
    {
      name: "Two trees",
      command: "tangle:1-1000 --no-axis --tanglegram before.nwk --against after.nwk",
      files:
        "=== before.nwk ===\n" +
        "((a:1,b:1):1,(c:1,d:1):1);\n" +
        "\n" +
        "=== after.nwk ===\n" +
        "((a:1,c:1):1,(b:1,d:1):1);\n",
    },
    {
      name: "Gene neighbourhoods",
      command: "ESX-1:1-4,000 --loci loci.bed --links hits.tsv --label 'ESX-1'",
      files:
        "=== loci.bed ===\n" +
        "H37Rv\t0\t1200\tespA\t0\t+\n" +
        "H37Rv\t1300\t2100\tespC\t0\t+\n" +
        "H37Rv\t2200\t3000\tespD\t0\t-\n" +
        "CDC1551\t0\t1200\tespA2\t0\t+\n" +
        "CDC1551\t2200\t3000\tespD2\t0\t-\n" +
        "\n" +
        "=== hits.tsv ===\n" +
        "espA\tespA2\t99.1\n" +
        "espD\tespD2\t97.4\n",
    },
    {
      name: "Protein domains",
      command: "protein:1-700 --domains domains.tsv --analysis Pfam",
      files:
        "=== domains.tsv ===\n" +
        "PknB\tmd5\t626\tPfam\tPF00069\tProtein kinase domain\t11\t275\t1e-40\tT\t01-01-2026\n" +
        "PknB\tmd5\t626\tPfam\tPF03793\tPASTA domain\t341\t400\t3e-10\tT\t01-01-2026\n" +
        "PknB\tmd5\t626\tPfam\tPF03793\tPASTA domain\t410\t468\t4e-10\tT\t01-01-2026\n" +
        "PknD\tmd5\t664\tPfam\tPF00069\tProtein kinase domain\t14\t277\t9e-40\tT\t01-01-2026\n" +
        "PknE\tmd5\t565\tPfam\tPF00069\tProtein kinase domain\t16\t280\t3e-39\tT\t01-01-2026\n",
    },
  ];

  var el = {};

  // A command line, not a shell. Quoting is handled because a title has spaces
  // in it and every example in the documentation shows one; everything else a
  // shell does is a shell's business and is not reimplemented here.
  function words(text) {
    var out = [];
    var word = "";
    var quote = null;
    var started = false;
    for (var i = 0; i < text.length; i++) {
      var c = text[i];
      if (quote) {
        if (c === quote) quote = null;
        else { word += c; started = true; }
      } else if (c === '"' || c === "'") {
        quote = c;
        started = true;
      } else if (c === "\\" && text[i + 1] === "\n") {
        i++;
      } else if (/\s/.test(c)) {
        if (started) { out.push(word); word = ""; started = false; }
      } else {
        word += c;
        started = true;
      }
    }
    if (started) out.push(word);
    return out;
  }

  // `=== name ===` on a line of its own starts a file. A header rather than a
  // separator, so a file may be empty and still be a file.
  function split(text) {
    var files = [];
    var current = null;
    text.split("\n").forEach(function (line) {
      var header = /^\s*={3,}\s*(.+?)\s*={3,}\s*$/.exec(line);
      if (header) {
        current = { name: header[1], body: [] };
        files.push(current);
      } else if (current) {
        current.body.push(line);
      }
    });
    return files.map(function (file) {
      return { name: file.name, body: file.body.join("\n") };
    });
  }

  function pack(argv, files) {
    var parts = [];
    var total = 0;
    function u32(n) {
      var b = new Uint8Array(4);
      new DataView(b.buffer).setUint32(0, n, true);
      parts.push(b);
      total += 4;
    }
    function str(s) {
      var bytes = encoder.encode(s);
      u32(bytes.length);
      parts.push(bytes);
      total += bytes.length;
    }
    u32(argv.length);
    argv.forEach(str);
    u32(files.length);
    files.forEach(function (file) { str(file.name); str(file.body); });

    var out = new Uint8Array(total);
    var at = 0;
    parts.forEach(function (part) { out.set(part, at); at += part.length; });
    return out;
  }

  function call(argv, files) {
    var input = pack(argv, files);
    var memory = wasm.memory;
    var inPtr = wasm.alloc(input.length);
    new Uint8Array(memory.buffer, inPtr, input.length).set(input);

    var outPtr = wasm.render(inPtr, input.length);
    wasm.dealloc(inPtr, input.length);

    // `memory.buffer` is detached and replaced whenever wasm grows its heap,
    // so it is read again here rather than reused from above.
    var view = new DataView(wasm.memory.buffer);
    var ok = view.getUint8(outPtr) === 1;
    var len = view.getUint32(outPtr + 1, true);
    var body = decoder.decode(new Uint8Array(wasm.memory.buffer, outPtr + 5, len));
    wasm.dealloc(outPtr, 5 + len);
    return { ok: ok, body: body };
  }

  function draw() {
    if (!wasm) return;
    var argv = words(el.command.value);
    var files = split(el.files.value);
    var answer;
    try {
      answer = call(argv, files);
    } catch (error) {
      answer = { ok: false, body: "the program stopped: " + error };
    }

    if (answer.ok) {
      el.output.className = "pg-output";
      el.output.innerHTML = answer.body;
      el.status.textContent = files.length + (files.length === 1 ? " file" : " files") + ", drawn here";
    } else {
      el.output.className = "pg-output pg-failed";
      el.output.textContent = "karyon: " + answer.body;
      el.status.textContent = "the command did not draw";
    }
  }

  function load(index) {
    var example = EXAMPLES[index];
    el.command.value = example.command;
    el.files.value = example.files;
    draw();
  }

  function start() {
    ["command", "files", "draw", "status", "output"].forEach(function (name) {
      el[name] = document.getElementById("pg-" + name);
    });
    var app = document.getElementById("karyon-playground-app");
    if (!app || !el.command) return;
    app.hidden = false;

    EXAMPLES.forEach(function (example, index) {
      var button = document.createElement("button");
      button.type = "button";
      button.className = "pg-chip";
      button.textContent = example.name;
      button.addEventListener("click", function () { load(index); });
      document.getElementById("pg-examples-buttons").appendChild(button);
    });

    el.draw.addEventListener("click", draw);
    // Ctrl or Cmd with Enter, which is what every editor on the page already
    // means by "run this".
    [el.command, el.files].forEach(function (box) {
      box.addEventListener("keydown", function (event) {
        if ((event.ctrlKey || event.metaKey) && event.key === "Enter") {
          event.preventDefault();
          draw();
        }
      });
    });

    el.command.value = EXAMPLES[0].command;
    el.files.value = EXAMPLES[0].files;
    el.draw.disabled = true;

    fetch(WASM)
      .then(function (response) {
        if (!response.ok) throw new Error(response.status + " fetching the program");
        return response.arrayBuffer();
      })
      // `instantiate` rather than `instantiateStreaming`, which needs the file
      // to arrive as application/wasm and fails outright when a host serves it
      // as something else.
      .then(function (bytes) { return WebAssembly.instantiate(bytes, {}); })
      .then(function (built) {
        wasm = built.instance.exports;
        el.draw.disabled = false;
        el.status.textContent = "ready";
        draw();
      })
      .catch(function (error) {
        el.status.textContent = "the program did not load: " + error.message;
        el.output.className = "pg-output pg-failed";
        el.output.textContent =
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
