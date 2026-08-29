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

  // The reference the sequence, ORF and pileup examples share, so the reads
  // carry the bases they are aligned to rather than a different string.
  function reference(bases) {
    var unit = "ACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAGG";
    var out = "";
    for (var i = 0; i < (bases || 300); i++) out += unit[i % 61];
    return out;
  }

  function fastaOf(name, seq) {
    var out = ">" + name + "\n";
    for (var i = 0; i < seq.length; i += 60) out += seq.slice(i, i + 60) + "\n";
    return out;
  }

  // A deterministic alignment with something in it, shared by the three
  // examples that read one. It was four sequences of sixty-one bases, pasted
  // out three times, and each of those figures was a picture of a track with
  // nothing to draw.
  //
  // What is in it and why. Three lineages, so a metadata strip has something to
  // say and a variable site can be diagnostic of one rather than scattered. A
  // conserved core where the letters of a logo stand full height and a
  // hypervariable stretch where they collapse, because a logo whose columns are
  // all one height is not showing anything. A deletion carried by one lineage,
  // so the alignment has gaps rather than being a rectangle of letters. And
  // sites that vary within a lineage as well as between them, so the
  // variable-sites panel is not just the lineage tree written out again.
  //
  // No Math.random. The same knobs give the same figure here and on the cards,
  // which is the promise the renderer itself makes.
  // One seeded generator for every example that needs a stream of numbers, so
  // the same knobs give the same figure here, on the cards and in a
  // screenshot. Math.random would make each of those a different picture.
  function rolls(seed) {
    var state = seed >>> 0;
    return function () {
      state = (state * 1103515245 + 12345) & 0x7fffffff;
      return state / 0x7fffffff;
    };
  }

  function alignment(rows, columns) {
    rows = rows || 18;
    columns = columns || 900;
    var bases = "ACGT";
    var seed = 20260829;
    function next() {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    }
    var core = [Math.round(columns * 0.34), Math.round(columns * 0.55)];
    var wild = [Math.round(columns * 0.66), Math.round(columns * 0.83)];
    var gap = [Math.round(columns * 0.13), Math.round(columns * 0.20)];

    var root = "";
    for (var c = 0; c < columns; c++) root += bases[Math.floor(next() * 4)];

    // Sites that tell one lineage from the others, spread over the whole length
    // so that panning finds them rather than one screen holding all of them.
    var marks = [[], [], []];
    for (var m = 0; m < 36; m++) marks[m % 3].push(Math.floor(next() * columns));

    var out = "";
    var names = [];
    for (var r = 0; r < rows; r++) {
      var lineage = r % 3;
      var seq = root.split("");
      for (var k = 0; k < marks[lineage].length; k++) {
        var at = marks[lineage][k];
        seq[at] = bases[(bases.indexOf(seq[at]) + 1) % 4];
      }
      for (var c2 = 0; c2 < columns; c2++) {
        var chance = c2 >= core[0] && c2 < core[1] ? 0.002
                   : c2 >= wild[0] && c2 < wild[1] ? 0.22
                   : 0.02;
        if (next() < chance) seq[c2] = bases[Math.floor(next() * 4)];
      }
      if (lineage === 1) {
        for (var g = gap[0]; g < gap[1]; g++) seq[g] = "-";
      }
      var name = "L" + (lineage + 1) + "_" + ("00" + (r + 1)).slice(-3);
      names.push(name);
      out += fastaOf(name, seq.join(""));
    }
    return { fasta: out, names: names, columns: columns };
  }

  // The sheet those names join to, drawn as strips beside the rows. What is
  // known about a sample is not a coordinate, so it is not a track of its own;
  // it hangs off the one that has the rows.
  function sheetFor(names, columns, lineageOf) {
    var out = "name\t" + columns.join("\t") + "\n";
    for (var i = 0; i < names.length; i++) {
      // The join is by name and never by row order, so the sheet has to carry
      // the names the track already knows. Where the name does not say which
      // group a row is in, the caller says.
      var lineage = lineageOf ? lineageOf(i, names[i]) : names[i].slice(0, 2);
      var cells = [];
      for (var c = 0; c < columns.length; c++) {
        var column = columns[c];
        cells.push(column === "lineage" ? lineage
                 : column === "source" ? (i % 3 === 0 ? "sputum" : i % 3 === 1 ? "culture" : "survey")
                 : column === "drug" ? (i % 4 === 0 ? "resistant" : "sensitive")
                 : column === "year" ? String(2014 + (i % 9))
                 : "");
      }
      out += names[i] + "\t" + cells.join("\t") + "\n";
    }
    return out;
  }

  // One alignment, written once and read twice: as ribbons between two
  // sequences and as a dot plot of the same blocks. The nine rows it used to
  // hold were all forward and all in order, so neither figure had anything in
  // it that the other could disagree about.
  //
  // What is in it now, in the order a reader meets it: a colinear opening, an
  // inversion, a stretch with no alignment at all, a translocation that lands
  // back near the start of the target, a resumption, and a duplication that
  // puts one piece of query on two places of target. Those are the four things
  // a synteny figure is for, and the dot plot shows the same four as a rising
  // diagonal, an anti-diagonal crossing it, a break, and a second diagonal
  // parallel to the first.
  function synteny(shape) {
    var next = rolls(400913);
    var rows = [];
    function block(qFrom, qTo, tFrom, tTo, strand, identity) {
      var span = Math.abs(qTo - qFrom);
      var matches = Math.round(span * identity);
      rows.push(["ctg1", 240000, qFrom, qTo, strand, "chrA", 300000,
                 Math.min(tFrom, tTo), Math.max(tFrom, tTo),
                 matches, span, 60].join("\t"));
    }
    var q = 0, t = 20000;
    for (var i = 0; i < 10; i++) {
      block(q, q + 5200, t, t + 5200, "+", 0.99 - i * 0.003);
      q += 6000; t += 6000;
    }
    if (shape !== "colinear only") {
      // An inversion: the query runs on while the target runs backwards.
      var it = 100000;
      for (var v = 0; v < 4; v++) {
        block(62000 + v * 6000, 62000 + v * 6000 + 5200, it, it - 5200, "-", 0.97);
        it -= 6000;
      }
      // Then nothing at all from 86,000 to 100,000 of the query.
      var tq = 100000, tt = 4000;
      for (var w = 0; w < 5; w++) {
        block(tq, tq + 5200, tt, tt + 5200, "+", 0.98);
        tq += 6000; tt += 6000;
      }
      var rq = 132000, rt = 152000;
      for (var r = 0; r < 6; r++) {
        block(rq, rq + 5200, rt, rt + 5200, "+", 0.96);
        rq += 6000; rt += 6000;
      }
    }
    if (shape === "everything") {
      // One piece of query against two places of target, which is what a
      // duplication looks like from the query's side.
      for (var d = 0; d < 3; d++) {
        block(170000 + d * 6000, 170000 + d * 6000 + 5200, 232000 + d * 6000, 232000 + d * 6000 + 5200, "+", 0.95);
        block(170000 + d * 6000, 170000 + d * 6000 + 5200, 60000 + d * 6000, 60000 + d * 6000 + 5200, "+", 0.93);
      }
      for (var x = 0; x < 6; x++) {
        var qq = 200000 + x * 5000;
        var tt2 = Math.floor(next() * 280000);
        block(qq, qq + 1400, tt2, tt2 + 1400, next() > 0.5 ? "+" : "-", 0.90);
      }
    }
    return rows.join("\n") + "\n";
  }

  var EXAMPLES = [
    {
      name: "A locus",
      bounds: { from: 756999, to: 766999, min: 270 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--aggregate", after: "--coverage", label: "How a pixel that covers many bases chooses", options: ["max", "mean", "min"] },
        { kind: "choice", flag: "--style", after: "--coverage", label: "The shape of the signal", options: ["area", "line", "bars"] },
        { kind: "toggle", flag: "--log", after: "--coverage", label: "A log scale for the depth" },
        { kind: "choice", flag: "--height", after: "--coverage",
          label: "How tall the depth track is", options: ["60", "110", "180"] },
      ],
      group: "Signal and annotation",
      command:
        "NC_000962.3:761,000-762,999 --coverage depth.bg --label depth --aggregate min \\\n" +
        "  --features genes.gff3 --label annotation \\\n" +
        "  --variants calls.vcf --label variants \\\n" +
        "  --title 'rpoB locus, resistance determining region'",
      // Wide enough that the window can be moved. The three rows that fitted
      // the region exactly refused on the first drag, which is the program
      // being right and the example being too small to show anything else.
      files: [
        { name: "depth.bg", body: "" },
        { name: "genes.gff3", body: "" },
        { name: "calls.vcf", body: "" },
      ],
      // Ten kilobases with something to see everywhere in them, so the window
      // can be taken anywhere inside its own bounds and still draw. Written
      // here rather than pasted because a reader wants to move it, not read it.
      make: function () {
        var depth = "", genes = "##gff-version 3\n", calls = "";
        var next = rolls(761000);
        // Per base rather than in hundred base bins. A bin wider than a pixel
        // does the aggregating before the program gets to, so --aggregate had
        // nothing left to decide; the deletion at 761,890 is 140 bases, which
        // is under a pixel at this width, and whether it shows at all is
        // exactly the question the flag answers.
        for (var at = 757000; at < 767000; at += 1) {
          var level = 52 + Math.round((next() - 0.5) * 14);
          if (at >= 761890 && at < 762030) level = 0;
          else if (at >= 759100 && at < 759400) level = Math.round(level * 0.45);
          else if (at >= 764200 && at < 764260) level = Math.round(level * 6.5);
          depth += "NC_000962.3 " + at + " " + (at + 1) + " " + level + "\n";
        }
        genes += "NC_000962.3 . gene 759807 763325 . + . Name=rpoB\n";
        genes += "NC_000962.3 . gene 761082 761162 . + . Name=RRDR\n";
        var others = [[757200, 758900, "rpoC", "+"], [763600, 765100, "rpsL", "-"],
                      [765400, 766800, "rrs", "+"], [758950, 759700, "Rv0666", "-"]];
        for (var g = 0; g < others.length; g++) {
          genes += "NC_000962.3 . gene " + others[g][0] + " " + others[g][1] +
                   " . " + others[g][3] + " . Name=" + others[g][2] + "\n";
        }
        var known = { 761109: 0.98, 761139: 0.55, 761155: 1.0, 761156: 0.21 };
        for (var v = 757060; v < 767000; v += 137) {
          var af = known[v] !== undefined ? known[v] : (0.05 + ((v % 17) / 20));
          var kind = af > 0.5 ? "missense_variant|MODERATE" : "synonymous_variant|LOW";
          calls += "NC_000962.3 " + v + " . C T . . AF=" + af.toFixed(2) +
                   ";ANN=T|" + kind + "|rpoB\n";
        }
        for (var k in known) {
          calls += "NC_000962.3 " + k + " . G A . . AF=" + known[k].toFixed(2) +
                   ";ANN=A|missense_variant|MODERATE|rpoB\n";
        }
        return [
          { name: "depth.bg", body: depth },
          { name: "genes.gff3", body: genes },
          { name: "calls.vcf", body: calls },
        ];
      },
    },
    {
      name: "A whole chromosome",
      bounds: { from: 1, to: 2000000, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--aggregate", after: "--coverage",
          label: "How a pixel that covers many bases chooses",
          options: ["max", "mean", "min"],
          says: "At this width one pixel covers two kilobases" },
        { kind: "toggle", flag: "--log", after: "--coverage", label: "A log scale for the depth" },
        { kind: "choice", flag: "--style", after: "--coverage",
          label: "The shape of the depth signal", options: ["area", "line", "bars"] },
      ],
      group: "Signal and annotation",
      // The depth used to be a sine wave with one dip in it, which made the
      // aggregate control nearly meaningless: a smooth curve looks much the
      // same whichever of the three a pixel picks. It is not smooth now. There
      // is a homozygous deletion at zero, a heterozygous loss at half depth, a
      // duplication, and eight single-bin repeat spikes at ten times baseline,
      // and those spikes are the whole point of the control: at two thousand
      // bases a pixel `max` draws them and `min` erases them, and neither is
      // wrong, which is why the program will not choose for you.
      command:
        "chr1:1-2,000,000 --coverage depth.bg --label depth --aggregate min \\\n" +
        "  --windows gc.bg --label 'GC content' --style steps",
      files: [
        { name: "depth.bg", body: "" },
        { name: "gc.bg", body: "" },
      ],
      make: function () {
        var next = rolls(120449);
        var depth = "", gc = "";
        // A hundred base bins, not a thousand. At nine hundred pixels a bin of
        // a thousand puts two of them under each pixel, and two numbers is not
        // enough for max, mean and min to disagree by much: the control was
        // there and it barely moved the picture. At a hundred bases a pixel
        // holds twenty-odd bins and the three answers come apart.
        var bins = 20000, wide = 100;
        var spikes = [1800, 4020, 6550, 8300, 10950, 13100, 15800, 18880];
        for (var i = 0; i < bins; i++) {
          var at = i * wide;
          var level = 38 + Math.round((next() - 0.5) * 18);
          if (i >= 9000 && i < 9600) level = 0;
          else if (i >= 12000 && i < 13200) level = Math.round(level * 0.5);
          else if (i >= 15000 && i < 15900) level = Math.round(level * 2.1);
          if (spikes.indexOf(i) >= 0) level = 430;
          depth += "chr1\t" + at + "\t" + (at + wide) + "\t" + level + "\n";
          if (i % 10 === 0) {
            gc += "chr1\t" + at + "\t" + (at + wide * 10) + "\t" +
                  (0.42 + 0.11 * Math.sin(i / 170) + (next() - 0.5) * 0.03).toFixed(3) + "\n";
          }
        }
        return [
          { name: "depth.bg", body: depth },
          { name: "gc.bg", body: gc },
        ];
      },
    },
    {
      name: "Two trees",
      bounds: { from: 1, to: 1000, min: 60 },
      controls: [
        { kind: "note", label: "Two phylogenies face to face, with a tie from each tip to the tip of the same name opposite. The ties are the figure; the trees are there to put the tips in the order each of them argues for. The axis is not a coordinate, so panning it would mean nothing." },
        { kind: "data", param: "tips", value: 20, label: "Taxa in both trees", options: [8, 20, 44] },
        { kind: "data", param: "disagreement", value: "both, and one taxon missing", label: "How much the two trees disagree",
          options: ["none", "a local swap", "a long jump", "both, and one taxon missing"] },
        { kind: "toggle", flag: "--no-axis", label: "The coordinate ruler, which measures nothing here" },
      ],
      group: "Phylogeny",
      // Four taxa cannot cross. The badge over the figure counts the crossings,
      // the linked tips and the unmatched ones, and with two tips a side there
      // was nothing for any of those three numbers to say.
      //
      // The two kinds of disagreement are drawn apart on purpose: a swap
      // between neighbouring tips is a short tie that crosses one other, and a
      // taxon that jumps across the figure is one long tie through everything,
      // and a tanglegram is read by telling those apart. The names are the join,
      // so a taxon only one tree has is drawn unmatched rather than dropped.
      command: "tangle:1-1000 --no-axis --tanglegram core.nwk --against accessory.nwk",
      files: [
        { name: "core.nwk", body: "" },
        { name: "accessory.nwk", body: "" },
      ],
      make: function (p) {
        var tips = (p && p.tips) || 20;
        var how = (p && p.disagreement) || "both, and one taxon missing";
        var names = [];
        for (var i = 0; i < tips; i++) {
          names.push("L" + (Math.floor(i / Math.ceil(tips / 4)) + 1) + "_" + ((i % Math.ceil(tips / 4)) + 1));
        }
        function newick(order, extra) {
          var groups = [];
          var per = Math.ceil(order.length / 4);
          for (var g = 0; g < order.length; g += per) {
            var members = order.slice(g, g + per).map(function (n) { return n + ":" + (0.4 + (n.length % 3) * 0.2).toFixed(2); });
            groups.push("(" + members.join(",") + "):1.0");
          }
          if (extra) groups.push(extra + ":1.6");
          return "(" + groups.join(",") + ");\n";
        }
        var left = names.slice();
        var right = names.slice();
        if (how === "a local swap" || how.indexOf("both") === 0) {
          var a = Math.floor(tips / 4), b = a + 1;
          var keep = right[a]; right[a] = right[b]; right[b] = keep;
        }
        if (how === "a long jump" || how.indexOf("both") === 0) {
          var from = Math.floor(tips / 2), to = tips - 2;
          var moved = right.splice(from, 1)[0];
          right.splice(to, 0, moved);
        }
        var only = how.indexOf("missing") >= 0 ? "L1_x" : null;
        return [
          { name: "core.nwk", body: newick(left, only) },
          { name: "accessory.nwk", body: newick(right, null) },
        ];
      },
    },
    {
      name: "Gene neighbourhoods",
      bounds: { from: 1, to: 19000, min: 1200 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "genomes", value: 5, label: "Genomes stacked",
          options: [2, 5, 9] },
        { kind: "data", param: "losses", value: "a deletion and an insertion each", label: "What the genomes are missing",
          options: ["nothing", "one deletion", "a deletion and an insertion each"] },
        { kind: "choice", flag: "--columns", after: "--traits",
          label: "Which strips to draw, in this order",
          options: ["clade", "clade,vaccine", "vaccine,year"] },
      ],
      group: "Comparisons across genomes",
      // Three rows of three genes could not show what a locus track is for. The
      // figure is read by looking down a column of orthologues and finding the
      // one row where it is not there, so the region has to hold enough genes
      // for a column to exist and enough rows for the gap to be a gap.
      //
      // This is the ESX-1 locus across five genomes. BCG has lost the six genes
      // of RD1, which is the deletion the vaccine strain is defined by, and the
      // strip beside the rows says which of them is a vaccine strain, so the
      // reader can check that the gap and the label agree.
      command: "ESX-1:1-9,000 --loci loci.bed --links hits.tsv --label 'ESX-1' --traits genomes.tsv --columns clade",
      files: [
        { name: "loci.bed", body: "" },
        { name: "hits.tsv", body: "" },
        { name: "genomes.tsv", body: "" },
      ],
      make: function (p) {
        var count = (p && p.genomes) || 5;
        var losses = (p && p.losses) || "a deletion and an insertion each";
        var genes = ["PE35", "PPE68", "esxB", "esxA", "espI", "eccD1", "espJ", "espK",
                     "espB", "mycP1", "eccB1", "eccCa1", "eccCb1", "eccE1", "espL", "espR"];
        var rd1 = ["PPE68", "esxB", "esxA", "espI", "eccD1", "espJ"];
        // BCG second on purpose. The losses control has nothing to say unless
        // a genome that lost something is on the page, and at two genomes it
        // was two complete rows and a knob that did nothing.
        var all = ["H37Rv", "BCG_Pasteur", "CDC1551", "Erdman", "M_microti",
                   "H37Ra", "Beijing", "Haarlem", "M_bovis"];
        var names = all.slice(0, count);
        var next = rolls(190822);
        var bed = "", links = "";
        for (var g = 0; g < names.length; g++) {
          var at = 600;
          var previous = null;
          for (var i = 0; i < genes.length; i++) {
            var gene = genes[i];
            var gone = losses !== "nothing" &&
                       ((names[g] === "BCG_Pasteur" && rd1.indexOf(gene) >= 0) ||
                        (names[g] === "M_microti" && ["esxB", "esxA", "espI"].indexOf(gene) >= 0) ||
                        (names[g] === "M_bovis" && gene === "espK"));
            var wide = 500 + Math.floor(next() * 700);
            if (!gone) {
              var here = gene + (g === 0 ? "" : "_" + (g + 1));
              bed += names[g] + "\t" + at + "\t" + (at + wide) + "\t" + here +
                     "\t0\t" + (i % 3 === 2 ? "-" : "+") + "\n";
              if (previous !== null) {
                links += previous + "\t" + here + "\t" + (94 + next() * 5.5).toFixed(1) + "\n";
              }
              previous = here;
            }
            at += wide + 200;
          }
          if (losses === "a deletion and an insertion each" && g > 0) {
            // Something only this genome has, so the unmatched outline is
            // exercised in the middle of the stack and not only at the bottom.
            bed += names[g] + "\t" + (7900 + g * 40) + "\t" + (9300 + g * 40) +
                   "\tIS6110_" + (g + 1) + "\t0\t+\n";
          }
        }
        // The links join gene to gene down the stack, so they are written from
        // one row to the next rather than from the first to all.
        var chain = "";
        for (var i2 = 0; i2 < genes.length; i2++) {
          for (var g2 = 0; g2 + 1 < names.length; g2++) {
            var from = genes[i2] + (g2 === 0 ? "" : "_" + (g2 + 1));
            var to = genes[i2] + "_" + (g2 + 2);
            if (bed.indexOf("\t" + from + "\t") < 0 || bed.indexOf("\t" + to + "\t") < 0) continue;
            chain += from + "\t" + to + "\t" + (94 + next() * 5.5).toFixed(1) + "\n";
          }
        }
        return [
          { name: "loci.bed", body: bed },
          { name: "hits.tsv", body: chain },
          { name: "genomes.tsv", body: "name\tclade\tvaccine\tyear\n" +
              names.map(function (n, i) {
                return n + "\t" + (i % 2 ? "animal" : "human") + "\t" +
                       (n.indexOf("BCG") === 0 ? "yes" : "no") + "\t" + (1998 + i * 3) + "\n";
              }).join("") },
        ];
      },
    },
    {
      name: "Protein domains",
      bounds: { from: 1, to: 1400, min: 100 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--analysis", after: "--domains",
          label: "Which member database to believe",
          options: ["Pfam", "SUPERFAMILY", "SMART", "PANTHER"] },
        { kind: "data", param: "proteins", value: 11, label: "Proteins stacked",
          options: [3, 11] },
        { kind: "choice", flag: "--columns", after: "--traits",
          label: "Which strips to draw, in this order",
          options: ["family", "family,essential", "essential"] },
      ],
      group: "Comparisons across genomes",
      // One protein and one database. The flag that picks a database had one
      // option, which is a control that cannot be wrong and therefore cannot
      // teach anything, and the shared residue axis had nothing to be shared
      // between.
      //
      // Eleven kinases across four databases is the figure the flag exists for:
      // the same protein gets different boundaries from different callers, and
      // switching between them is the fastest way to see that a domain call is
      // an opinion rather than a coordinate.
      command: "protein:1-700 --domains domains.tsv --analysis Pfam --traits proteins.tsv --columns family",
      files: [
        { name: "domains.tsv", body: "" },
        { name: "proteins.tsv", body: "" },
      ],
      make: function (p) {
        var count = (p && p.proteins) || 11;
        var kinases = [
          ["PknA", 431], ["PknB", 626], ["PknD", 664], ["PknE", 565], ["PknF", 472],
          ["PknG", 750], ["PknH", 626], ["PknI", 585], ["PknJ", 583], ["PknK", 1317],
          ["PknL", 429],
        ].slice(0, count);
        var next = rolls(551903);
        var out = "";
        function row(name, length, database, accession, label, from, to) {
          out += [name, "md5", length, database, accession, label, from, to,
                  "1e-" + (10 + Math.floor(next() * 40)), "T", "01-01-2026"].join("\t") + "\n";
        }
        for (var k = 0; k < kinases.length; k++) {
          var name = kinases[k][0], length = kinases[k][1];
          // The same kinase domain, called four ways, with the boundaries
          // deliberately disagreeing by a few dozen residues.
          row(name, length, "Pfam", "PF00069", "Protein kinase domain", 11, 275);
          row(name, length, "SUPERFAMILY", "SSF56112", "Protein kinase-like", 4, 291);
          row(name, length, "SMART", "SM00220", "S_TKc", 17, 268);
          row(name, length, "PANTHER", "PTHR24347", "Serine/threonine-protein kinase", 1, 320);
          if (name === "PknB" || name === "PknE" || name === "PknJ") {
            var repeats = name === "PknB" ? 4 : 1;
            for (var r = 0; r < repeats; r++) {
              row(name, length, "Pfam", "PF03793", "PASTA domain", 341 + r * 70, 400 + r * 70);
              row(name, length, "SMART", "SM00740", "PASTA", 336 + r * 70, 405 + r * 70);
            }
          }
          if (name === "PknG") {
            row(name, length, "Pfam", "PF00301", "Rubredoxin", 74, 140);
            row(name, length, "Pfam", "PF13424", "TPR repeat", 480, 660);
            row(name, length, "PANTHER", "PTHR24347", "Kinase, TPR region", 470, 700);
          }
          if (name === "PknK") {
            row(name, length, "Pfam", "PF13191", "AAA domain", 420, 620);
            row(name, length, "Pfam", "PF00196", "LuxR-type HTH", 1180, 1270);
            row(name, length, "SUPERFAMILY", "SSF52540", "P-loop NTPase", 405, 640);
          }
          if (name === "PknD") row(name, length, "Pfam", "PF13360", "PQQ-like propeller", 400, 640);
          if (name === "PknH") row(name, length, "Pfam", "PF14589", "Extracellular domain", 420, 560);
        }
        return [
          { name: "domains.tsv", body: out },
          { name: "proteins.tsv", body: "name\tfamily\tessential\n" +
              kinases.map(function (k, i) {
                return k[0] + "\t" + (i % 3 === 0 ? "PASTA" : i % 3 === 1 ? "soluble" : "membrane") +
                       "\t" + (k[0] === "PknA" || k[0] === "PknB" ? "yes" : "no") + "\n";
              }).join("") },
        ];
      },
    },
    {
      name: "One molecule at a time",
      bounds: { from: 900, to: 4900, min: 400 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--context", after: "--bisulfite",
          label: "Which cytosine context was called",
          options: ["CpG", "CHG", "CHH"] },
        { kind: "data", param: "molecules", value: 24, label: "Molecules sequenced",
          options: [8, 24, 60] },
        { kind: "choice", flag: "--max-rows", after: "--bisulfite",
          label: "How many rows before it stops and counts the rest",
          options: ["8", "40", "all"] },
        { kind: "toggle", flag: "--no-names", after: "--bisulfite",
          label: "Leave out the name beside each row" },
        { kind: "toggle", flag: "--no-axis", label: "Leave out the coordinate ruler" },
      ],
      group: "Reads and molecules",
      // Ten reads over nine sites, all of them called in one context, is a
      // picture with nothing in it: half the molecules were methylated
      // everywhere and half nowhere, and --context had one setting because the
      // file held one context.
      //
      // What is here now is an imprinting control region: an island where the
      // molecules split cleanly into methylated and unmethylated stripes, and a
      // shore outside it where the same molecules are confetti. That contrast
      // is what the track exists to show, and it needs a window wide enough to
      // hold both, which is why the region runs to 4,900 rather than 200.
      //
      // The molecules start and end in different places, so the left and right
      // edges of the grid are ragged. An open circle is a site the molecule
      // covered and did not modify, and no circle at all is a site it never
      // reached, and those are two different statements that a rectangular grid
      // cannot tell apart.
      command: "chr11:1,400-3,400 --bisulfite calls.txt --context CpG --label 'H19 ICR'",
      files: [{ name: "calls.txt", body: "" }],
      make: function (p) {
        var molecules = (p && p.molecules) || 24;
        var next = rolls(3400291);
        var island = [1005, 2600];
        var sites = [];
        for (var at = 1005; at < 4800; ) {
          sites.push(at);
          at += at < island[1] ? 18 + Math.floor(next() * 24) : 60 + Math.floor(next() * 80);
        }
        var chg = [], chh = [];
        for (var c = 950; c < 4880; c += 150) chg.push(c);
        for (var h = 930; h < 4880; h += 90) chh.push(h);

        var out = "";
        for (var m = 0; m < molecules; m++) {
          var name = "read" + (m + 1);
          var methylated = m % 2 === 0;
          var from = 950 + Math.floor(next() * 700);
          var to = 4200 + Math.floor(next() * 650);
          function row(pos, on) {
            // Column two is the case of the call written out again, and the
            // reader refuses the file where the two disagree, so they are
            // written from one decision rather than two.
            out += name + "\t" + (on ? "+" : "-") + "\tchr11\t" + pos + "\t";
          }
          for (var i = 0; i < sites.length; i++) {
            var pos = sites[i];
            if (pos < from || pos > to) continue;
            var inside = pos >= island[0] && pos < island[1];
            var on = inside ? (methylated ? next() > 0.05 : next() < 0.05)
                            : next() < 0.35;
            row(pos, on);
            out += (on ? "Z" : "z") + "\n";
          }
          for (var g = 0; g < chg.length; g++) {
            if (chg[g] < from || chg[g] > to) continue;
            var gon = next() < 0.04;
            row(chg[g], gon);
            out += (gon ? "X" : "x") + "\n";
          }
          for (var e = 0; e < chh.length; e++) {
            if (chh[e] < from || chh[e] > to) continue;
            var eon = next() < 0.02;
            row(chh[e], eon);
            out += (eon ? "H" : "h") + "\n";
          }
        }
        return [{ name: "calls.txt", body: out }];
      },
    },
    {
      name: "Reference bases",
      bounds: { from: 1, to: 1500, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "orf", value: "one long reading frame", label: "What is in the sequence",
          options: ["one long reading frame", "stops everywhere", "two frames overlapping"] },
        { kind: "toggle", flag: "--no-region-label",
          label: "Leave out the locus printed at the top right" },
      ],
      group: "Signal and annotation",
      // Three hundred bases of a repeating unit gave six frame lanes that were
      // all the same, because a sequence that repeats every sixty-one bases has
      // stops in the same places in every frame. There is a designed open
      // reading frame in here now, so one lane carries a long bar and the other
      // five are broken up, which is the difference the track is drawn to show.
      //
      // The window matters too: below about five pixels a base the letters give
      // way to coloured blocks and then to a single line, which is the sequence
      // track choosing what it can afford rather than drawing letters nobody
      // could read.
      command: "chr1:150-750 --axis --sequence ref.fa --label reference --orfs ref.fa --label 'reading frames'",
      files: [{ name: "ref.fa", body: "" }],
      make: function (p) {
        var telling = (p && p.orf) || "one long reading frame";
        var next = rolls(140277);
        var bases = "ACGT";
        var stops = ["TAA", "TAG", "TGA"];
        function filler(n, stopRate) {
          var out = "";
          for (var i = 0; i < n; i += 3) {
            if (next() < stopRate) out += stops[Math.floor(next() * 3)];
            else {
              var codon = "";
              for (var b = 0; b < 3; b++) codon += bases[Math.floor(next() * 4)];
              out += stops.indexOf(codon) >= 0 ? "CTG" : codon;
            }
          }
          return out.slice(0, n);
        }
        function orf(codons) {
          var out = "ATG";
          for (var c = 1; c < codons - 1; c++) {
            var codon = "";
            for (var b = 0; b < 3; b++) codon += bases[Math.floor(next() * 4)];
            out += stops.indexOf(codon) >= 0 ? "GCT" : codon;
          }
          return out + "TAA";
        }
        var seq;
        if (telling === "stops everywhere") {
          seq = filler(1500, 0.5);
        } else if (telling === "two frames overlapping") {
          seq = filler(299, 0.35) + orf(120) + "A" + orf(90) + filler(1500, 0.35);
        } else {
          seq = filler(299, 0.35) + orf(134) + filler(1500, 0.35);
        }
        return [{ name: "ref.fa", body: fastaOf("chr1", seq.slice(0, 1500)) }];
      },
    },
    {
      name: "A cytogenetic ideogram",
      bounds: { from: 1, to: 46709983, min: 100000 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--height", after: "--ideogram",
          label: "How thick the chromosome is drawn", options: ["18", "34", "70"] },
        { kind: "toggle", flag: "--no-region-label",
          label: "Leave out the locus printed at the top right" },
      ],
      group: "Signal and annotation",
      // Ten bands of equal width over two megabases drew a ramp, which is what
      // a cytoBand table looks like when nobody gave it a chromosome. This is
      // the shape of a real acrocentric: a short arm of stalk and variable
      // heterochromatin, an off-centre waist, and a long arm whose bands run
      // from under two megabases to over eight. Panning it is the point, since
      // a band is the one annotation whose width is the whole of its meaning.
      command: "chr21:1-46,709,983 --ideogram bands.txt --label chr21",
      files: [{ name: "bands.txt", body: "" }],
      make: function () {
        var bands = [
          [0, 3100000, "p13", "gvar"], [3100000, 6800000, "p12", "stalk"],
          [6800000, 10900000, "p11.2", "gvar"], [10900000, 11700000, "p11.1", "acen"],
          [11700000, 12400000, "q11.1", "acen"], [12400000, 14300000, "q11.2", "gneg"],
          [14300000, 20200000, "q21.1", "gpos100"], [20200000, 24100000, "q21.2", "gneg"],
          [24100000, 28800000, "q21.3", "gpos75"], [28800000, 32400000, "q22.11", "gneg"],
          [32400000, 36100000, "q22.12", "gpos50"], [36100000, 38400000, "q22.13", "gneg"],
          [38400000, 46709983 - 8300000, "q22.2", "gpos25"],
          [46709983 - 8300000, 46709983, "q22.3", "gneg"],
        ];
        var out = "";
        for (var i = 0; i < bands.length; i++) {
          out += "chr21\t" + bands[i][0] + "\t" + bands[i][1] + "\t" +
                 bands[i][2] + "\t" + bands[i][3] + "\n";
        }
        return [{ name: "bands.txt", body: out }];
      },
    },
    {
      name: "A genome-wide scan",
      bounds: { from: 1, to: 1000000, min: 5000 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "markers", value: 2500, label: "Markers tested",
          options: [500, 2500, 8000], says: "More tests, and the null climbs with them" },
        { kind: "data", param: "signal", value: "two loci and some singletons", label: "What is under the peaks",
          options: ["nothing but noise", "one locus", "two loci and some singletons"] },
        { kind: "choice", flag: "--threshold", after: "--manhattan",
          label: "The line it is read against",
          options: ["5", "genome-wide", "9"],
          says: "genome-wide is a correction for a million tests" },
        { kind: "choice", flag: "--height", after: "--manhattan",
          label: "How tall the scan is drawn", options: ["80", "160", "300"] },
      ],
      group: "Association and genotype",
      // The point of the marker control is the thing a Manhattan plot is
      // usually read wrongly for: with five hundred draws the highest point of
      // pure noise is around three, and with eight thousand it is over four, so
      // the height of the tallest bar means nothing until you know how many
      // tests were run. Turning the knob moves the noise floor while the two
      // real peaks stay where they are.
      command: "chr1:1-1,000,000 --manhattan assoc.tsv --label association",
      files: [{ name: "assoc.tsv", body: "" }],
      make: function (p) {
        var markers = (p && p.markers) || 2500;
        var signal = (p && p.signal) || "two loci and some singletons";
        var next = rolls(880021);
        var step = Math.floor(1000000 / markers);
        var out = "pos\tneglog10p\n";
        for (var i = 0; i < markers; i++) {
          var at = 1 + i * step;
          var u = Math.max(1e-5, next());
          var value = -Math.log(u) / Math.LN10;
          if (signal !== "nothing but noise") {
            var d1 = (at - 615000) / 6000;
            value = Math.max(value, 11.5 * Math.exp(-d1 * d1 / 2));
          }
          if (signal === "two loci and some singletons") {
            var d2 = (at - 240000) / 2500;
            value = Math.max(value, 8.2 * Math.exp(-d2 * d2 / 2));
            if (i % Math.max(1, Math.floor(markers / 4)) === 7) value = 5.6 + next() * 0.8;
          }
          out += at + "\t" + value.toFixed(3) + "\n";
        }
        return [{ name: "assoc.tsv", body: out }];
      },
    },
    {
      name: "A genotype matrix",
      bounds: { from: 1, to: 30500, min: 800 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "samples", value: 32, label: "Samples in the matrix",
          options: [12, 32, 80] },
        { kind: "data", param: "structure", value: "three clades", label: "What the samples are",
          options: ["one population", "three clades", "three clades and a mixture"] },
        { kind: "choice", flag: "--columns", after: "--traits",
          label: "Which strips to draw, in this order",
          options: ["lineage", "lineage,drug", "drug,year,source", "source"] },
      ],
      group: "Association and genotype",
      // Twelve rows by seven columns of alternating 0.95 and 0.05 was a
      // chequerboard, and a chequerboard is the one pattern a genotype matrix
      // cannot mean anything by. The sites are irregular now, the way variant
      // sites are, and the samples fall into blocks that a strip beside them
      // names, so the question the figure answers is whether the blocks in the
      // matrix agree with the labels on the rows.
      command: "chr1:1-12,000 --matrix geno.tsv --label 'allele fraction' --traits samples.tsv --columns lineage",
      files: [
        { name: "geno.tsv", body: "" },
        { name: "samples.tsv", body: "" },
      ],
      make: function (p) {
        var samples = (p && p.samples) || 32;
        var structure = (p && p.structure) || "three clades";
        var next = rolls(720113);
        var sites = [];
        for (var at = 150; at < 30400; at += 200 + Math.floor(next() * 140)) sites.push(at);

        var names = [];
        for (var i = 0; i < samples; i++) names.push("ERR" + ("00" + (i + 1)).slice(-3));
        var clade = function (i) {
          if (structure === "one population") return 0;
          var per = Math.ceil(samples / 3);
          return Math.floor(i / per);
        };
        var out = "sample\t" + sites.join("\t") + "\n";
        for (var r = 0; r < samples; r++) {
          var mine = clade(r);
          var mixed = structure.indexOf("mixture") >= 0 && r % 11 === 5;
          var row = [names[r]];
          for (var c = 0; c < sites.length; c++) {
            var derived = (c % 3) === mine || (mine === 1 && c >= 28 && c < 48);
            if (mixed) derived = next() < 0.5;
            var value = derived ? 0.9 + next() * 0.1 : next() * 0.1;
            row.push(next() < 0.02 ? "NA" : value.toFixed(2));
          }
          out += row.join("\t") + "\n";
        }
        return [
          { name: "geno.tsv", body: out },
          { name: "samples.tsv", body: sheetFor(names, ["lineage", "drug", "year", "source"],
              function (i) { return "L" + (clade(i) + 1); }) },
        ];
      },
    },
    {
      name: "A read pileup",
      bounds: { from: 1, to: 1600, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "copies", value: 8, label: "Reads starting at each step",
          options: [1, 3, 8, 20], unit: "x",
          says: "The track stacks 40 rows and counts the rest" },
        { kind: "data", param: "reach", value: "a mix", label: "How far each read reaches",
          options: ["55 bases", "120 bases", "a mix"] },
        { kind: "data", param: "carrying", value: "one variant",
          label: "What the reads disagree about",
          options: ["nothing", "one variant", "a strand artefact"],
          says: "A base is coloured only where it differs from the reference" },
        { kind: "choice", flag: "--max-rows", after: "--pileup",
          label: "How deep it is drawn before it counts the rest",
          options: ["10", "40", "120", "all"] },
        { kind: "data", param: "quality", value: "a mix", label: "Mapping quality",
          options: ["all confident", "a mix"],
          says: "Only a mix gives --fade-by-mapq anything to fade" },
        { kind: "toggle", flag: "--fade-by-mapq", after: "--pileup",
          label: "Draw a read fainter the less sure the aligner was" },
        { kind: "note", label: "The fade is on to begin with, because a pileup where every read is equally solid is a pileup that has not been asked the question. Turn it off and the two mapping qualities draw the same figure, which is the point of the knob above it." },
        { kind: "choice", flag: "--row-height", after: "--pileup",
          label: "How tall one read is", options: ["6", "12", "24"] },
      ],
      group: "Reads and molecules",
      command:
        "chr1:1-400 --sequence ref.fa --label reference \\\n" +
        "  --pileup reads.sam --with-sequence ref.fa --fade-by-mapq --label reads",
      files: [
        { name: "ref.fa", body: "" },
        { name: "reads.sam", body: "" },
      ],
      // What a pileup answers to is how many reads there are and how far each
      // one reaches, so those are the two controls. Coverage deep enough to
      // reach the forty row cap is the thing worth watching: the track stops
      // opening rows there and counts the rest, rather than drawing off the
      // bottom of the figure.
      //
      // The reference is given twice on purpose: once as a track of its own,
      // so the letters are on the page, and once to the pileup, which is what
      // lets it colour a base that disagrees. Until --with-sequence reached a
      // pileup the second of those was impossible and every read agreed, which
      // is a picture of the track doing nothing.
      //
      // A variant every read carries is a variant. One only the forward reads
      // carry is the shape of an artefact, and telling those two apart is what
      // drawing the reads is for.
      make: function (p) {
        var ref = reference(1600);
        var copies = p && p.copies ? p.copies : 3;
        var reach = (p && p.reach) || "55 bases";
        var carrying = (p && p.carrying) || "one variant";
        var mixedQuality = (p && p.quality) !== "all confident";
        var swap = { A: "G", C: "T", G: "A", T: "C" };
        var sites = carrying === "nothing" ? [] : [97, 184, 320];
        var forwardOnly = carrying === "a strand artefact";
        var sam = "";
        var n = 0;
        // Which reads carry the variant is drawn rather than counted. Keying
        // it on the read number put it in step with the strand whenever the
        // depth was even, since both alternate, and then "one variant" and "a
        // strand artefact" wrote the same file and drew the same figure.
        var next = rolls(902117);
        for (var at = 1; at < 1540; at += 6) {
          for (var copy = 0; copy < copies; copy++) {
            var len = reach === "120 bases" ? 120
                    : reach === "a mix" ? 30 + ((at + copy * 13) % 110)
                    : 55 + ((at + copy * 7) % 12);
            if (at - 1 + len > ref.length) len = ref.length - at + 1;
            if (len < 5) continue;
            var forward = (copy % 2) === 0;
            var seq = ref.slice(at - 1, at - 1 + len).split("");
            for (var v = 0; v < sites.length; v++) {
              var off = sites[v] - at;
              if (off < 0 || off >= len) continue;
              var carries = forwardOnly ? forward : next() < 0.5;
              if (carries) seq[off] = swap[seq[off]] || "N";
            }
            // A read the aligner could have placed anywhere carries a low
            // mapping quality, and the whole point of --fade-by-mapq is to
            // stop it looking as solid as one that could not. With every read
            // at 60 the flag has nothing to do, so the knob above says whether
            // they vary.
            var mapq = mixedQuality ? [60, 2, 30, 0][n % 4] : 60;
            sam += "r" + (n++) + "\t" + (forward ? 0 : 16) + "\tchr1\t" + at +
                   "\t" + mapq + "\t" + len + "M\t*\t0\t0\t" + seq.join("") + "\t*\n";
          }
        }
        return [
          { name: "ref.fa", body: fastaOf("chr1", ref) },
          { name: "reads.sam", body: sam },
        ];
      },
    },
    {
      name: "One molecule in pieces",
      bounds: { from: 1, to: 12000, min: 400 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "molecules", value: 44, label: "Molecules in the file",
          options: [12, 44, 90] },
        { kind: "data", param: "kinds", value: "everything, including three-piece molecules", label: "What the pieces did",
          options: ["forward hops only", "hops and inversions", "everything, including three-piece molecules"] },
        { kind: "choice", flag: "--height", after: "--split-reads",
          label: "How much room the rows get", options: ["80", "160", "300"] },
      ],
      group: "Reads and molecules",
      // Thirty molecules that all hopped forwards, all in two pieces, all the
      // same distance. A split read is worth drawing because the pieces can go
      // backwards, land on the other strand, or come in threes, and none of
      // those was in the file.
      command: "chr1:1-5,000 --split-reads split.sam --label 'split reads'",
      files: [{ name: "split.sam", body: "" }],
      make: function (p) {
        var molecules = (p && p.molecules) || 44;
        var kinds = (p && p.kinds) || "everything, including three-piece molecules";
        var next = rolls(311017);
        var sam = "";
        for (var i = 0; i < molecules; i++) {
          var a = 300 + Math.floor(i * (11200 / molecules));
          var which = kinds === "forward hops only" ? 0
                    : kinds === "hops and inversions" ? (i % 2)
                    : (i % 4);
          if (which === 2) {
            // Three pieces, so "segment 2 of 3" exists at all.
            var b = a + 1400 + Math.floor(next() * 900);
            var c = b + 1300 + Math.floor(next() * 700);
            sam += "m" + i + "\t0\tchr1\t" + a + "\t60\t400M900S\t*\t0\t0\t*\t*\tSA:Z:chr1," +
                   b + ",+,400S450M450S,60,0;chr1," + c + ",+,850S450M,60,0;\n";
            continue;
          }
          var back = which === 3;
          var far = back ? a - 1500 - Math.floor(next() * 600) : a + 1600 + Math.floor(next() * 900);
          if (far < 1) far = a + 1600;
          var reverse = which === 1;
          sam += "m" + i + "\t0\tchr1\t" + a + "\t60\t600M700S\t*\t0\t0\t*\t*\t" +
                 "SA:Z:chr1," + far + "," + (reverse ? "-" : "+") + ",600S700M,60,0;\n";
        }
        return [{ name: "split.sam", body: sam }];
      },
    },
    {
      name: "Modified bases",
      bounds: { from: 1, to: 10000, min: 300 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--modification", after: "--methylation",
          label: "Which modified base was counted", options: ["m", "h"] },
        { kind: "choice", flag: "--min-reads", after: "--methylation",
          label: "The fewest reads behind a site for it to be drawn",
          options: ["5", "20", "45"] },
        { kind: "choice", flag: "--height", after: "--methylation",
          label: "How much room the two strand lanes get",
          options: ["60", "120", "220"] },
        { kind: "data", param: "depth", value: "mixed", label: "How deeply the sites were covered",
          options: ["thin", "mixed", "deep"],
          says: "Calls under 5x are dropped and counted" },
        { kind: "toggle", flag: "--no-axis", label: "Leave out the coordinate ruler" },
      ],
      group: "Reads and molecules",
      // The file used to hold one modification code, one strand and one
      // coverage, so the lower lane was empty, the strand pair was never seen,
      // the fade never varied and the note the track prints when it drops thin
      // calls never appeared. Four of the things the track does could not be
      // seen at all.
      //
      // Now both strands are called, the shore between 4,000 and 4,800 is
      // hemimethylated so the two lanes disagree, and 5hmC is a different
      // picture rather than the same one recoloured. The coverage knob is what
      // makes the 5x filter visible: on "thin" the figure says how many calls
      // it left out.
      command: "chr1:1,500-5,000 --methylation calls.bed --modification m --label '5mC'",
      files: [{ name: "calls.bed", body: "" }],
      make: function (p) {
        var depth = (p && p.depth) || "mixed";
        var pool = depth === "thin" ? [2, 3, 4, 5, 7, 9]
                 : depth === "deep" ? [28, 40, 55, 70, 90, 120]
                 : [3, 4, 6, 12, 28, 40, 55, 70];
        var next = rolls(778211);
        var out = "";
        for (var at = 60; at < 9900; at += 30 + Math.floor(next() * 60)) {
          var island = at >= 2000 && at < 4000;
          var shore = at >= 4000 && at < 4800;
          for (var strand = 0; strand < 2; strand++) {
            var reverse = strand === 1;
            for (var code = 0; code < 2; code++) {
              var which = code === 0 ? "m" : "h";
              var base = which === "h" ? 8
                       : island ? 88
                       : shore ? (reverse ? 10 : 55)
                       : 9;
              var pct = Math.max(0, Math.min(100, base + (next() - 0.5) * 12));
              var cov = pool[Math.floor(next() * pool.length)];
              out += "chr1\t" + at + "\t" + (at + 1) + "\t" + which + "\t" + cov + "\t" +
                     (reverse ? "-" : "+") + "\t" + at + "\t" + (at + 1) + "\t0,0,0\t" +
                     cov + "\t" + pct.toFixed(2) + "\t1\t1\t0\t0\t0\t0\t0\n";
            }
          }
        }
        return [{ name: "calls.bed", body: out }];
      },
    },
    {
      name: "Alignment ribbons",
      bounds: { from: 1, to: 240000, min: 8000 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "shape", value: "everything", label: "What the two sequences did",
          options: ["colinear only", "with an inversion and a jump", "everything"] },
        { kind: "choice", flag: "--height", after: "--synteny",
          label: "How much room the ribbons get", options: ["80", "160", "300"] },
      ],
      group: "Alignments and rearrangements",
      // Nine forward blocks in order was a picture of two sequences agreeing,
      // which is the one case where ribbons have nothing to say. The alignment
      // is the same one the dot plot below reads, so the two figures can be
      // opened side by side and compared, which is the reason both are here.
      command: "ctg1:1-190,000 --synteny aln.paf --label 'against chrA'",
      files: [{ name: "aln.paf", body: "" }],
      make: function (p) {
        return [{ name: "aln.paf", body: synteny((p && p.shape) || "everything") }];
      },
    },
    {
      name: "The same PAF as a dot plot",
      bounds: { from: 1, to: 240000, min: 8000 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "shape", value: "everything", label: "What the two sequences did",
          options: ["colinear only", "with an inversion and a jump", "everything"] },
      ],
      group: "Alignments and rearrangements",
      // The same file as the ribbons above, which is the point of having both.
      // A rising diagonal is colinearity; the anti-diagonal crossing it is the
      // inversion; the diagonal that restarts near the bottom is the piece that
      // moved; and two diagonals at the same query is one piece of it living in
      // two places of the target.
      command: "ctg1:1-240,000 --dotplot aln.paf --label 'dot plot'",
      files: [{ name: "aln.paf", body: "" }],
      make: function (p) {
        return [{ name: "aln.paf", body: synteny((p && p.shape) || "everything") }];
      },
    },
    {
      name: "A multiple alignment",
      bounds: { from: 1, to: 900, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "rows", value: 18, label: "Sequences in the alignment",
          options: [6, 18, 40, 90], says: "Every one is a row, and nothing caps them" },
        { kind: "data", param: "columns", value: 900, label: "Columns it is long",
          options: [300, 900, 2400] },
        { kind: "choice", flag: "--style", after: "--msa",
          label: "Only the cells that disagree, or all of them",
          options: ["differences", "all"] },
        { kind: "choice", flag: "--compare-to", after: "--msa",
          label: "The row the others are read against",
          options: ["L1_001", "L2_002", "L3_003"],
          says: "The consensus without it" },
        { kind: "choice", flag: "--max-rows", after: "--msa",
          label: "How many rows before it stops and counts the rest",
          options: ["8", "40", "all"] },
        { kind: "toggle", flag: "--no-names", after: "--msa",
          label: "Leave out the name beside each row" },
        { kind: "choice", flag: "--columns", after: "--traits",
          label: "Which strips to draw, in this order",
          options: ["lineage", "lineage,drug", "drug,year,source", "source"] },
      ],
      group: "Sequence alignment",
      // The window is the thing to move. Zoomed out the alignment is a block of
      // agreement with a few columns of disagreement running down it, and the
      // deletion the second lineage carries is a hole in the middle of it.
      // Zoomed in far enough the letters arrive, which is the track deciding
      // for itself what it can afford to draw.
      command: "aln:1-400 --msa aln.fa --label alignment --traits meta.tsv --columns lineage",
      files: [
        { name: "aln.fa", body: "" },
        { name: "meta.tsv", body: "" },
      ],
      make: function (p) {
        var made = alignment(p && p.rows, p && p.columns);
        return [
          { name: "aln.fa", body: made.fasta },
          { name: "meta.tsv", body: sheetFor(made.names, ["lineage", "drug", "year", "source"]) },
        ];
      },
    },
    {
      name: "Only the variable sites",
      bounds: { from: 1, to: 900, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "rows", value: 18, label: "Sequences in the alignment",
          options: [6, 18, 40, 90] },
        { kind: "data", param: "columns", value: 300, label: "Columns the alignment is long",
          options: [150, 300, 600], says: "The panel keeps only the ones that disagree" },
        { kind: "choice", flag: "--compare-to", after: "--snps",
          label: "The row the others are read against",
          options: ["L1_001", "L2_002", "L3_003"],
          says: "Whichever record came first, without it" },
        { kind: "toggle", flag: "--no-counts", after: "--snps",
          label: "Leave out the difference count down the right" },
        { kind: "choice", flag: "--max-rows", after: "--snps",
          label: "How many rows before it stops and counts the rest",
          options: ["8", "40", "all"] },
        { kind: "toggle", flag: "--no-names", after: "--snps",
          label: "Leave out the name beside each row" },
        { kind: "choice", flag: "--columns", after: "--traits",
          label: "Which strips to draw, in this order",
          options: ["lineage", "lineage,drug", "drug,year,source", "source"] },
      ],
      group: "Sequence alignment",
      // The same alignment as the one above, which is why both are here. This
      // panel throws away every column where nothing disagrees, so hundreds of
      // columns become a few dozen and the three lineages fall into stripes
      // that the strips beside them name.
      command: "aln:1-400 --snps aln.fa --label 'variable sites' --traits meta.tsv --columns lineage",
      files: [
        { name: "aln.fa", body: "" },
        { name: "meta.tsv", body: "" },
      ],
      make: function (p) {
        var made = alignment(p && p.rows, p && p.columns);
        return [
          { name: "aln.fa", body: made.fasta },
          { name: "meta.tsv", body: sheetFor(made.names, ["lineage", "drug", "year", "source"]) },
        ];
      },
    },
    {
      name: "A sequence logo",
      bounds: { from: 1, to: 900, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "rows", value: 18, label: "Sequences the logo is built from",
          options: [6, 18, 40, 90],
          says: "Few sequences make every column look certain" },
        { kind: "data", param: "columns", value: 900, label: "Columns it is long",
          options: [300, 900, 2400] },
      ],
      group: "Sequence alignment",
      // The alignment has a conserved core and a hypervariable stretch on
      // purpose, because that difference is the whole of what a logo says. Move
      // the window from one to the other and the letters go from full height to
      // a scatter of stubs. The sequence count is the other half of it: six
      // sequences agreeing looks exactly like certainty, and ninety agreeing
      // is certainty, and the logo is drawn from a count that cannot tell them
      // apart unless somebody says how many there were.
      command: "aln:1-400 --logo aln.fa --label 'information content'",
      files: [{ name: "aln.fa", body: "" }],
      make: function (p) {
        var made = alignment(p && p.rows, p && p.columns);
        return [{ name: "aln.fa", body: made.fasta }];
      },
    },
    {
      name: "A phylogeny",
      bounds: { from: 1, to: 100, min: 60 },
      controls: [
        { kind: "note", label: "A phylogeny has no coordinates, so there is nothing to pan across. What it has is a shape, and support values that live in the tooltips: hover a node." },
        { kind: "data", param: "tips", value: 24, label: "Taxa in the tree", options: [8, 24, 60] },
        { kind: "data", param: "support", value: "bootstraps", label: "Node support",
          options: ["none", "bootstraps"], says: "Hover an internal node to read it" },
        { kind: "toggle", flag: "--no-axis", label: "The coordinate ruler, which measures nothing here" },
      ],
      group: "Phylogeny",
      // Six tips and no support values. The tree layout puts support in the
      // tooltips rather than on the branches by default, which is the right
      // choice for a printed figure and means the reader of a live page has
      // something to hover, so long as the file carries any. This one does.
      command: "tree:1-100 --no-axis --tree tree.nwk --label phylogeny",
      files: [{ name: "tree.nwk", body: "" }],
      make: function (p) {
        var tips = (p && p.tips) || 24;
        var carrying = (p && p.support) !== "none";
        var next = rolls(66041);
        var supports = [98, 100, 72, 89, 55, 100, 94, 63, 100, 81, 77, 96];
        var used = 0;
        function label() { return carrying ? String(supports[used++ % supports.length]) : ""; }
        var per = Math.ceil(tips / 4);
        var groups = [];
        for (var g = 0; g < 4; g++) {
          var members = [];
          for (var i = 0; i < per && g * per + i < tips; i++) {
            var name = "L" + (g + 1) + "_" + ("0" + (i + 1)).slice(-2);
            members.push(name + ":" + (0.2 + next() * 0.5).toFixed(3));
          }
          // Nest half of each group one level deeper, so the tree has a shape
          // rather than being a comb.
          var half = Math.max(1, Math.floor(members.length / 2));
          var inner = "(" + members.slice(0, half).join(",") + ")" + label() + ":0.4";
          var rest = members.slice(half);
          groups.push("(" + [inner].concat(rest).join(",") + ")" + label() + ":0.8");
        }
        return [{ name: "tree.nwk", body: "(" + groups.join(",") + ")" + label() + ";\n" }];
      },
    },
    {
      name: "Recombination on a tree",
      bounds: { from: 1, to: 60000, min: 2000 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "blocks", value: 21, label: "Blocks the caller found",
          options: [6, 21, 48] },
        { kind: "data", param: "clades", value: "one block that is not a clade", label: "Who carries them",
          options: ["clades only", "one block that is not a clade"],
          says: "A block whose taxa are not a clade is drawn row by row" },
        { kind: "choice", flag: "--columns", after: "--traits",
          label: "Which strips to draw, in this order",
          options: ["lineage", "lineage,source", "source"] },
      ],
      group: "Phylogeny",
      // Four blocks over six taxa. The whole claim a clade track makes is that
      // one rectangle stands for one event on one branch, and the module spends
      // most of itself refusing to make that claim when the taxa carrying a
      // block are not a clade. With four blocks there was never a case where it
      // had to refuse, so the reader never saw the thing it is careful about.
      command: "NC_011900.1:1-20,000 --clades gubbins.gff --with-tree tree.nwk --label recombination --traits taxa.tsv --columns lineage",
      files: [
        { name: "gubbins.gff", body: "" },
        { name: "tree.nwk", body: "" },
        { name: "taxa.tsv", body: "" },
      ],
      make: function (p) {
        var count = (p && p.blocks) || 21;
        var awkward = (p && p.clades) !== "clades only";
        var next = rolls(220517);
        var names = [];
        for (var g = 0; g < 4; g++) {
          for (var i = 0; i < 5; i++) names.push("L" + (g + 1) + "_" + ("0" + (i + 1)).slice(-2));
        }
        var groups = [];
        for (var q = 0; q < 4; q++) {
          var members = names.slice(q * 5, q * 5 + 5).map(function (n) {
            return n + ":" + (0.2 + next() * 0.4).toFixed(3);
          });
          groups.push("(" + members.join(",") + ")" + (70 + Math.floor(next() * 30)) + ":0.9");
        }
        var tree = "(" + groups.join(",") + ");\n";

        var gff = "##gff-version 3\n";
        var at = 800;
        for (var b = 0; b < count; b++) {
          var wide = 600 + Math.floor(next() * 2300);
          var lineage = b % 4;
          var size = 2 + Math.floor(next() * 4);
          var taxa = names.slice(lineage * 5, lineage * 5 + Math.min(5, size));
          if (awkward && b === 1) {
            // Two taxa from two different lineages, which is not a clade. The
            // track will not draw one rectangle over them.
            taxa = [names[0], names[12]];
          }
          gff += "SEQUENCE\tGUBBINS\tCDS\t" + at + "\t" + (at + wide) +
                 "\t0.000\t.\t0\tnode=\"N" + b + "\";taxa=\"" + taxa.join(" ") + "\";\n";
          at += wide + 200 + Math.floor(next() * 1200);
          if (at > 58000) at = 800 + Math.floor(next() * 400);
        }
        return [
          { name: "gubbins.gff", body: gff },
          { name: "tree.nwk", body: tree },
          { name: "taxa.tsv", body: sheetFor(names, ["lineage", "source"]) },
        ];
      },
    },
    {
      name: "Structural variants",
      bounds: { from: 1, to: 200000, min: 4000 },
      controls: [
        { kind: "region" },
        { kind: "data", param: "kinds", value: "everything, including a breakend", label: "What the caller found",
          options: ["deletions only", "the four span types", "everything, including a breakend"] },
        { kind: "data", param: "support", value: "a wide range", label: "Reads behind each call",
          options: ["thin, 4 to 12", "a wide range", "deep, 40 to 90"],
          says: "The arc's weight is the support, so a flat range draws flat arcs" },
        { kind: "choice", flag: "--style", after: "--coverage",
          label: "The shape of the depth profile", options: ["area", "line", "bars"] },
      ],
      group: "Alignments and rearrangements",
      // Forty calls of four types at even spacing, with support drawn from a
      // narrow range, drew forty arcs of much the same weight. The weight is
      // the measurement, so a narrow range is a figure with the measurement
      // taken out of it, and that is what the support control is for.
      //
      // The depth track underneath is the other half. A deletion should be a
      // hole in the coverage and a duplication a step up in it, and an arc that
      // does not line up with one is the reason to doubt the call.
      command:
        "chr1:1-70,000 --coverage depth.bg --label depth \\\n" +
        "  --structural sv.vcf --label 'structural variants'",
      files: [
        { name: "depth.bg", body: "" },
        { name: "sv.vcf", body: "" },
      ],
      make: function (p) {
        var kinds = (p && p.kinds) || "everything, including a breakend";
        var range = (p && p.support) || "a wide range";
        var next = rolls(660311);
        var calls = [
          ["DEL", 4000, 9200], ["DUP", 14000, 19500], ["INV", 24000, 31000],
          ["INS", 34000, 34850], ["DEL", 41000, 42200], ["BND", 52000, 141000],
          ["DUP", 61000, 88000], ["INV", 96000, 99000], ["DEL", 104000, 152000],
          ["INS", 120000, 120600], ["DEL", 160000, 161500], ["INV", 166000, 168500],
          ["DUP", 172000, 178000], ["DEL", 182000, 186400], ["DEL", 190000, 196000],
        ];
        var vcf = "";
        var copy = [];
        for (var i = 0; i < calls.length; i++) {
          var type = calls[i][0], from = calls[i][1], to = calls[i][2];
          if (kinds === "deletions only" && type !== "DEL") continue;
          if (kinds !== "everything, including a breakend" && type === "BND") continue;
          var reads = range === "thin, 4 to 12" ? 4 + Math.floor(next() * 9)
                    : range === "deep, 40 to 90" ? 40 + Math.floor(next() * 51)
                    : 5 + Math.floor(next() * 80);
          if (type === "BND") {
            vcf += "chr1\t" + from + "\tbnd" + i + "\tN\tt[chr1:" + to + "[\t.\t.\tSVTYPE=BND;SUPPORT=" + reads + "\n";
          } else {
            vcf += "chr1\t" + from + "\tsv" + i + "\tN\t<" + type + ">\t.\t.\tSVTYPE=" + type +
                   ";END=" + to + ";SVLEN=" + (to - from) + ";SUPPORT=" + reads + "\n";
          }
          if (type === "DEL" || type === "DUP") copy.push([from, to, type]);
        }
        // Depth that agrees with the calls, so an arc has something under it to
        // be read against.
        var depth = "";
        for (var b = 0; b < 200000; b += 500) {
          var level = 42;
          for (var c = 0; c < copy.length; c++) {
            if (b >= copy[c][0] && b < copy[c][1]) level = copy[c][2] === "DEL" ? 3 : 78;
          }
          depth += "chr1\t" + b + "\t" + (b + 500) + "\t" + (level + Math.floor(next() * 8)) + "\n";
        }
        return [
          { name: "depth.bg", body: depth },
          { name: "sv.vcf", body: vcf },
        ];
      },
    },
    {
      name: "Segmented copy number",
      bounds: { from: 1, to: 2000000, min: 20000 },
      controls: [
        { kind: "note", label: "A caller's segments, drawn against the ploidy the flag names. Balanced is where the ladder says it is, not where the data averages out, which is why --ploidy is required and not guessed." },
        { kind: "region" },
        { kind: "choice", flag: "--sample", after: "--copy-number",
          label: "Which sample of the table to draw",
          options: ["diagnosis", "relapse"] },
        { kind: "choice", flag: "--ploidy", after: "--copy-number",
          label: "Where balanced sits on the ladder", options: ["2", "3", "4"] },
      ],
      group: "Signal and annotation",
      // Seven segments in one sample could not show the two things this track
      // is for. There are two samples now, from the same patient at two times,
      // so the flag that picks one has something to pick between and the reader
      // can watch an amplification deepen from six copies to nine and a loss go
      // from one copy to none. A segment with cn1 at zero is drawn as
      // heterozygosity lost, which is a different statement from a plain loss
      // and needs both samples to be seen as a change.
      command: "chr8:1-2,000,000 --copy-number segments.cns --ploidy 2 --sample diagnosis --label 'copy number'",
      files: [{ name: "segments.cns", body: "" }],
      make: function () {
        var next = rolls(90124);
        var edges = [0];
        while (edges[edges.length - 1] < 2000000) {
          edges.push(Math.min(2000000, edges[edges.length - 1] + 40000 + Math.floor(next() * 180000)));
        }
        var out = "chromosome\tstart\tend\tgene\tlog2\tcn\tcn1\tcn2\tdepth\tweight\tsample\n";
        function write(sample, plan) {
          for (var i = 0; i + 1 < edges.length; i++) {
            var from = edges[i], to = edges[i + 1];
            var cn = plan(from, i);
            var na = i === 4;
            var log2 = na ? "NA" : (Math.log(Math.max(0.05, cn / 2)) / Math.LN2).toFixed(4);
            var minor = cn === 0 ? 0 : cn <= 1 ? 0 : Math.floor(cn / 2);
            out += "chr8\t" + from + "\t" + to + "\t" + (i === 6 ? "MYC" : "-") + "\t" +
                   log2 + "\t" + (na ? "NA" : cn) + "\t" + (na ? "NA" : minor) + "\t" +
                   (na ? "NA" : cn - minor) + "\t" +
                   (30 + Math.floor(next() * 40)) + "\t" + (0.4 + next() * 0.6).toFixed(2) +
                   "\t" + sample + "\n";
          }
        }
        write("diagnosis", function (from, i) {
          if (i >= 6 && i <= 8) return 6;
          if (i >= 11 && i <= 13) return 1;
          if (i >= 16) return 3;
          return 2;
        });
        write("relapse", function (from, i) {
          if (i >= 6 && i <= 8) return 9;
          if (i >= 11 && i <= 13) return 0;
          if (i >= 16) return 2;
          return 2;
        });
        return [{ name: "segments.cns", body: out }];
      },
    },
    {
      name: "Per-base model attribution",
      bounds: { from: 1, to: 4000, min: 60 },
      controls: [
        { kind: "note", label: "The bases themselves are the bars. A letter above the line is one the model leaned on, and one below is one it pulled away from, so the height is a signed score rather than a count." },
        { kind: "region" },
        { kind: "choice", flag: "--height", after: "--dynseq",
          label: "How tall the band is, which is how tall a letter can get",
          options: ["120", "200", "320"] },
        { kind: "data", param: "gaps", value: "two", label: "Stretches the model never scored",
          options: ["none", "two"], says: "A gap is a break in the rule, not a zero" },
      ],
      group: "Sequence alignment",
      // A hundred and twenty bases is one screen, so the window slider had
      // nowhere to go and the track's own ladder was never climbed. Over four
      // thousand bases it is: zoomed out the figure is an envelope, part way in
      // it is coloured bars, and close up the letters arrive, and those are
      // three regimes rather than three sizes of the same drawing.
      //
      // The gaps are the other half. A base the model was never run over is
      // absent from the file, and the track draws a break in the rule rather
      // than a zero, because a score of nothing and a score of zero are two
      // different statements about the model.
      command:
        "promoter:600-1,400 --dynseq attribution.bg --with-sequence promoter.fa \\\n" +
        "  --label attribution",
      files: [
        { name: "attribution.bg", body: "" },
        { name: "promoter.fa", body: "" },
      ],
      make: function (p) {
        var holes = (p && p.gaps) !== "none";
        var next = rolls(410903);
        var bases = "ACGT";
        var motifs = [
          [742, "TATAAAAG", 1.0], [1180, "CACGTG", 0.72], [1875, "GGGGCGGGGC", -0.85],
          [2460, "TGACTCA", 0.64], [3120, "CAGGTG", -0.55], [3540, "TATAAAAG", 0.88],
        ];
        var seq = "";
        for (var i = 0; i < 4000; i++) seq += bases[Math.floor(next() * 4)];
        var letters = seq.split("");
        for (var m = 0; m < motifs.length; m++) {
          var at = motifs[m][0], word = motifs[m][1];
          for (var c = 0; c < word.length; c++) letters[at + c] = word[c];
        }
        seq = letters.join("");

        var gaps = holes ? [[600, 660], [2800, 2960]] : [];
        var scores = "";
        for (var b = 0; b < 4000; b++) {
          var skip = false;
          for (var g = 0; g < gaps.length; g++) {
            if (b >= gaps[g][0] && b < gaps[g][1]) skip = true;
          }
          if (skip) continue;
          var value = 0.05 * Math.sin(b / 37) + 0.03 * Math.sin(b / 9.5);
          for (var k = 0; k < motifs.length; k++) {
            var from = motifs[k][0], width = motifs[k][1].length, peak = motifs[k][2];
            if (b >= from && b < from + width) {
              value = peak * (0.75 + 0.25 * Math.sin((b - from) / width * Math.PI));
            }
          }
          scores += "promoter\t" + b + "\t" + (b + 1) + "\t" + value.toFixed(4) + "\n";
        }
        return [
          { name: "attribution.bg", body: scores },
          { name: "promoter.fa", body: fastaOf("promoter", seq) },
        ];
      },
    },
    {
      name: "Splice junctions",
      bounds: { from: 1, to: 15500, min: 1500 },
      controls: [
        { kind: "note", label: "An arc per intron, thicker for the junctions more reads crossed. Multi-mapping reads are counted separately and never added in: a read that mapped in four places is one read." },
        { kind: "region" },
        { kind: "choice", flag: "--style", after: "--coverage",
          label: "The shape of the depth profile", options: ["area", "line", "bars"] },
        { kind: "toggle", flag: "--log", after: "--coverage",
          label: "A log scale, which is what shows the introns are not empty" },
        { kind: "choice", flag: "--min-reads", after: "--junctions",
          label: "The fewest reads across an intron for its arc to be drawn",
          options: ["50", "200", "500"] },
        { kind: "toggle", flag: "--no-counts", after: "--junctions",
          label: "Leave out the read count over each arc" },
        { kind: "choice", flag: "--height", after: "--junctions",
          label: "How much room the arcs get to miss each other",
          options: ["60", "120", "240"] },
        { kind: "data", param: "splicing", value: "the full picture", label: "What the gene is doing",
          options: ["one transcript", "an exon skipped", "the full picture"] },
      ],
      group: "Reads and molecules",
      // Five arcs over empty space was not a sashimi plot, it was a drawing of
      // arcs. The figure the field reads junctions in is arcs over the coverage
      // they came from, because the question is always which of two arcs over
      // one exon carries more reads, and an arc alone cannot be compared to
      // anything.
      //
      // So there are two tracks now, and the gene has eleven exons with enough
      // room between them that two arcs fit in view at once. The weight
      // encoding also needed the room: with five junctions the stroke widths
      // spanned 1.58x for a count range of 7.4x, which is a measurement the
      // reader cannot see.
      command:
        "chr1:1-9,500 --coverage depth.bg --label 'RNA-seq depth' \\\n" +
        "  --junctions SJ.out.tab --label junctions",
      files: [
        { name: "depth.bg", body: "" },
        { name: "SJ.out.tab", body: "" },
      ],
      make: function (p) {
        var telling = (p && p.splicing) || "the full picture";
        var next = rolls(51221);
        var exons = [];
        var at = 300;
        for (var e = 0; e < 11; e++) {
          var wide = 300 + Math.floor(next() * 400);
          exons.push([at, at + wide]);
          at += wide + 800 + Math.floor(next() * 500);
        }
        var depth = "";
        var last = exons[exons.length - 1][1];
        for (var b = 0; b < 15500; b += 25) {
          var over = 5;
          for (var x = 0; x < exons.length; x++) {
            if (b >= exons[x][0] && b < exons[x][1]) {
              over = 900 - x * 60 + Math.floor(next() * 60);
              break;
            }
          }
          depth += "chr1\t" + b + "\t" + (b + 25) + "\t" + over + "\n";
        }
        // STAR's columns: chrom, first intron base, last intron base, strand,
        // motif, annotated, unique reads, multi-mapping reads, overhang.
        var sj = "";
        function junction(from, to, unique, multi, annotated) {
          sj += "chr1\t" + (from + 1) + "\t" + to + "\t1\t2\t" + (annotated ? 1 : 0) +
                "\t" + unique + "\t" + multi + "\t" + (28 + Math.floor(next() * 12)) + "\n";
        }
        var counts = [910, 780, 620, 530, 455, 340, 300, 210, 160, 120];
        for (var j = 0; j < exons.length - 1; j++) {
          junction(exons[j][1], exons[j + 1][0], counts[j % counts.length], 0, true);
        }
        if (telling !== "one transcript") {
          // Exon three skipped, which is the arc that has to be read against
          // the two it spans rather than on its own.
          junction(exons[1][1], exons[3][0], 96, 0, true);
          junction(exons[8][1], exons[10][0], 143, 0, true);
        }
        if (telling === "the full picture") {
          // An alternative donor a hundred and forty bases inside exon five,
          // and one junction nobody has annotated, carried mostly by reads that
          // mapped in more than one place.
          junction(exons[4][0] + 140, exons[5][0], 41, 0, true);
          junction(exons[6][1] + 120, exons[7][0] - 60, 7, 61, false);
        }
        return [
          { name: "depth.bg", body: depth },
          { name: "SJ.out.tab", body: sj },
        ];
      },
    },
  ];

  // ---------------------------------------------------------------------
  // The controls, which belong to the example
  // ---------------------------------------------------------------------
  //
  // Not the same for every one, because the figures are not the same thing.
  // A window to slide is what a signal over a chromosome has and a tanglegram
  // has not, and offering the tanglegram one anyway would be a control that
  // does nothing and says nothing about why.
  //
  // No example offers the same flag twice, and that is a limit rather than a
  // choice. `after` says which track an option belongs to, but it only steers
  // where a new word is inserted: setFlag, flagOf and hasFlag in karyon-wasm.js
  // all find the flag with argv.indexOf and take the first one in the command.
  // Setting --style on a variants track therefore deletes the --style the
  // coverage track above it was carrying, and both selects read back the same
  // value. Two --style controls in one example need those three made aware of
  // `after` first.
  //
  // Most controls are a flag, so turning one rewrites a word of the command and
  // the command stays the thing that decides. The `data` kind is the exception,
  // and it is here because of what the command line can and cannot say. It has
  // 55 flags; the library behind it has some three hundred builder options, and
  // 242 of those have no flag at all. Outside coverage, windows and variants
  // there is almost nothing a flag can change, which is why twenty of these two
  // dozen examples had a window slider and nothing else beside it.
  //
  // What those tracks answer to is the file. So a `data` control hands a value
  // to the example's make(), which writes the file again, and the reader
  // watches a pileup fill up or an alignment gain rows rather than reading that
  // it would. It rewrites a file instead of a word of the command, and that is
  // not a sleight of hand: the file tabs are editable, so it does in one
  // gesture what the reader can already do by hand, and the command above goes
  // on saying exactly what drew the figure.

  var current = null;
  var params = {};

  function bounds() {
    return current && current.bounds ? current.bounds : null;
  }

  function retune() {
    // Written from the command rather than remembered, so a control and the
    // text above it cannot come apart when one of them is edited by hand.
    var strip = el.controls;
    if (!current || !current.controls) { strip.textContent = ""; return; }
    strip.textContent = "";

    current.controls.forEach(function (spec) {
      if (spec.kind === "note") {
        var note = document.createElement("p");
        note.className = "pg-note";
        note.textContent = spec.label;
        strip.appendChild(note);
        return;
      }

      var row = document.createElement("label");
      row.className = "pg-control";

      var name = document.createElement("span");
      name.className = "pg-control-name";
      strip.appendChild(row);

      if (spec.kind === "region") {
        var edge = bounds();
        var where = K.locus(el.command.value);
        if (!edge || !where) return;
        var span = where.end - where.start + 1;
        var room = edge.to - edge.from + 1;

        name.textContent = "Window";
        row.appendChild(name);

        var slide = document.createElement("input");
        slide.type = "range";
        slide.min = String(edge.from);
        slide.max = String(Math.max(edge.from, edge.to - span + 1));
        slide.value = String(Math.min(Math.max(where.start, edge.from), edge.to - span + 1));
        slide.step = String(Math.max(1, Math.round(span / 50)));
        slide.setAttribute("aria-label", "Where the window sits");
        slide.addEventListener("input", function () {
          var at = parseInt(slide.value, 10);
          move(K.within(K.retarget(el.command.value, at, at + span - 1), edge));
        });
        row.appendChild(slide);

        var wide = document.createElement("label");
        wide.className = "pg-control";
        var zname = document.createElement("span");
        zname.className = "pg-control-name";
        zname.textContent = "Width";
        wide.appendChild(zname);

        var zoom = document.createElement("input");
        zoom.type = "range";
        // A log scale, because a window runs from tens of bases to millions
        // and a linear slider spends nine tenths of itself on the last order
        // of magnitude.
        zoom.min = "0";
        zoom.max = "1000";
        var lo = Math.log(Math.max(K.MIN_SPAN, edge.min || 60, Math.min(room, 60)));
        var hi = Math.log(room);
        zoom.value = String(Math.round(((Math.log(span) - lo) / Math.max(1e-9, hi - lo)) * 1000));
        zoom.setAttribute("aria-label", "How many bases are in view");
        zoom.addEventListener("input", function () {
          var want = Math.round(Math.exp(lo + (parseInt(zoom.value, 10) / 1000) * (hi - lo)));
          var here = K.locus(el.command.value);
          var middle = here.start + (here.end - here.start) / 2;
          move(K.within(K.retarget(el.command.value, middle - want / 2, middle + want / 2 - 1), edge));
        });
        wide.appendChild(zoom);

        var says = document.createElement("output");
        says.className = "pg-control-says";
        says.textContent = K.grouped(span) + " bases";
        wide.appendChild(says);
        strip.appendChild(wide);
        return;
      }

      if (spec.kind === "data") {
        name.textContent = spec.label || spec.param;
        row.appendChild(name);

        var vary = document.createElement("select");
        spec.options.forEach(function (option) {
          var o = document.createElement("option");
          o.value = String(option);
          o.textContent = String(option) + (spec.unit || "");
          if (String(option) === String(params[spec.param])) o.selected = true;
          vary.appendChild(o);
        });
        vary.setAttribute("aria-label", spec.label || spec.param);
        vary.addEventListener("change", function () {
          var value = Number(vary.value);
          params[spec.param] = vary.value !== "" && !isNaN(value) ? value : vary.value;
          regenerate();
        });
        row.appendChild(vary);

        if (spec.says) {
          var told = document.createElement("span");
          told.className = "pg-control-says";
          told.textContent = spec.says;
          row.appendChild(told);
        }
        return;
      }

      name.textContent = spec.flag;
      name.title = spec.label || "";
      row.appendChild(name);

      if (spec.kind === "toggle") {
        var box = document.createElement("input");
        box.type = "checkbox";
        box.checked = K.hasFlag(el.command.value, spec.flag);
        box.addEventListener("change", function () {
          move(K.setFlag(el.command.value, spec.flag, box.checked ? true : null, spec.after));
        });
        row.appendChild(box);
      } else if (spec.kind === "choice") {
        var pick = document.createElement("select");
        var now = K.flagOf(el.command.value, spec.flag);
        // The empty option is what the flag not being there looks like, and
        // that is a real state: the program has a default and says so.
        [""].concat(spec.options).forEach(function (option) {
          var o = document.createElement("option");
          o.value = option;
          o.textContent = option || "(the default)";
          if (option === (now === null ? "" : now)) o.selected = true;
          pick.appendChild(o);
        });
        pick.addEventListener("change", function () {
          move(K.setFlag(el.command.value, spec.flag, pick.value || null, spec.after));
        });
        row.appendChild(pick);
      }
      if (spec.label) {
        var why = document.createElement("span");
        why.className = "pg-control-says";
        why.textContent = spec.label;
        row.appendChild(why);
      }
    });
  }

  // ---------------------------------------------------------------------
  // Drawing
  // ---------------------------------------------------------------------

  var drawn = null;

  function draw() {
    if (!K.ready()) return;
    save();
    var answer = K.run(el.command.value, files, el.plot.clientWidth - 24);

    if (answer.ok) {
      // Rebuilt rather than reassigned: writing the whole attribute took
      // `pg-dragging` off one pixel into a drag, so the grabbing cursor
      // flickered back to a hand for the rest of it.
      el.plot.classList.remove("pg-failed");
      el.plot.classList.add("pg-plot");
      el.plot.classList.toggle("pg-live", el.live.checked);
      // A figure is one image and is announced as one. A refusal is words, and
      // `role="img"` would make them presentational, so the pane only wears
      // that role while it holds a picture.
      el.plot.setAttribute("role", "img");
      el.plot.setAttribute("aria-label", "The figure this command draws");
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
      retune();
    } else {
      el.plot.className = "pg-plot pg-failed";
      el.plot.removeAttribute("role");
      el.plot.removeAttribute("aria-label");
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
    // Reachable only while it does something, so a reader tabbing past a
    // static figure is not stopped at it for nothing.
    el.plot.tabIndex = on ? 0 : -1;
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
    move(K.within(K.panned(el.command.value, moved / Math.max(1, el.plot.clientWidth)), bounds()));
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
    // Only once the figure has focus, so a wheel over a figure nobody has
    // clicked scrolls the page rather than being eaten by it.
    if (!el.live.checked || document.activeElement !== el.plot) return;
    event.preventDefault();
    var box = el.plot.getBoundingClientRect();
    var at = Math.min(1, Math.max(0, (event.clientX - box.left) / box.width));
    move(K.within(K.zoomed(el.command.value, event.deltaY > 0 ? 1.25 : 0.8, at), bounds()));
  }

  // ---------------------------------------------------------------------
  // Files, as tabs
  // ---------------------------------------------------------------------

  // The file the editor is holding, and the file itself rather than where it
  // sits. Two reasons. It is nothing at all until a file has been put in it,
  // and saving without that wrote an empty editor over the first file of every
  // example the moment it loaded. And a position is not an identity: closing
  // the open tab shifted every file after it down one, so the save that ran on
  // the way out wrote the closed file's text into whichever file had taken its
  // place. Closing `depth.bg` left `genes.gff3` holding a bedGraph.
  var showing = null;

  function save() {
    if (showing && files.indexOf(showing) >= 0) showing.body = el.file.value;
  }

  function show(index) {
    save();
    active = Math.max(0, Math.min(index, files.length - 1));
    el.file.value = files.length ? files[active].body : "";
    el.file.disabled = !files.length;
    showing = files.length ? files[active] : null;
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
      // The close control lives inside the tab, so without this the tab's own
      // name is read out as "depth.bg times".
      tab.setAttribute("aria-label", file.name);
      tab.title = file.name + "  (F2 renames, Delete removes)";
      tab.addEventListener("click", function () { show(index); });

      var rename = function () {
        var name = prompt("Name of this file, as the command calls it", file.name);
        if (name) { file.name = name.trim(); tabs(); draw(); }
      };
      var drop = function () {
        files.splice(index, 1);
        show(Math.min(active, files.length - 1));
        draw();
        var next = el.tabs.querySelector(".pg-tab");
        if (next) next.focus();
      };

      // A double click renames it, which is what a pointer expects, and F2
      // does the same, which is what every file list in the world uses: the
      // name is what the command calls the file by, and renaming it was bound
      // to a gesture a keyboard cannot make.
      tab.addEventListener("dblclick", rename);
      tab.addEventListener("keydown", function (event) {
        if (event.key === "F2") { event.preventDefault(); rename(); }
        else if (event.key === "Delete") { event.preventDefault(); drop(); }
      });

      // A button rather than a span, so it is in the tab order and answers to
      // Enter: a file could not be removed without a pointer at all.
      var shut = document.createElement("button");
      shut.type = "button";
      shut.className = "pg-shut";
      shut.textContent = "×";
      shut.setAttribute("aria-label", "Remove " + file.name);
      shut.addEventListener("click", function (event) {
        event.stopPropagation();
        drop();
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

  // Writes the example's files again from the values its `data` controls hold,
  // leaving the command alone. The reader may have panned somewhere, and losing
  // that on every turn of a knob would make the knob not worth turning.
  function regenerate() {
    if (!current || !current.make) return;
    showing = null;
    files = current.make(params).map(function (file) {
      return { name: file.name, body: file.body };
    });
    show(Math.min(active, files.length - 1));
    draw();
    retune();
  }

  function load(example) {
    current = example;
    params = {};
    (example.controls || []).forEach(function (spec) {
      if (spec.kind !== "data") return;
      params[spec.param] = spec.value !== undefined ? spec.value : spec.options[0];
    });
    el.command.value = example.command;
    // The editor is holding the last example's file, not this one's, so it
    // has nothing to save.
    showing = null;
    files = (example.make ? example.make(params) : example.files).map(function (file) {
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

  // ---------------------------------------------------------------------
  // The example picker
  // ---------------------------------------------------------------------

  var opened = false;
  var lastFocus = null;

  var previews = 0;

  /* Every figure the program draws names its own title and description with the
     same two ids, which is right for a document holding one figure and wrong
     for this panel, which holds two dozen. Injected as they come, all of them
     point their `aria-labelledby` at the first card's title, so a screen reader
     is told the same thing about all of them. Each one gets its own prefix on
     the way in, on every id and on every reference to one. */
  function unique(root, prefix) {
    var nodes = root.querySelectorAll("[id]");
    var i;
    for (i = 0; i < nodes.length; i++) {
      nodes[i].id = prefix + nodes[i].id;
    }
    var all = root.querySelectorAll("*");
    for (i = 0; i < all.length; i++) {
      var el = all[i];
      var labels = el.getAttribute("aria-labelledby");
      if (labels) {
        el.setAttribute(
          "aria-labelledby",
          labels.split(/\s+/).map(function (id) { return prefix + id; }).join(" ")
        );
      }
      for (var a = 0; a < el.attributes.length; a++) {
        var attr = el.attributes[a];
        if (attr.value.indexOf("url(#") >= 0) {
          attr.value = attr.value.replace(/url\(#/g, "url(#" + prefix);
        } else if (attr.name === "href" || attr.name === "xlink:href") {
          if (attr.value.charAt(0) === "#") attr.value = "#" + prefix + attr.value.slice(1);
        }
      }
    }
  }

  // What an example opens as, which is what its card should show. Reading it
  // from the controls rather than letting make() fall back to its own defaults
  // is the only thing that keeps the card and the figure the same picture.
  function opensAs(example) {
    var chosen = {};
    (example.controls || []).forEach(function (spec) {
      if (spec.kind !== "data") return;
      chosen[spec.param] = spec.value !== undefined ? spec.value : spec.options[0];
    });
    return chosen;
  }

  function preview(example, into) {
    // Drawn by the program, out of the same files the example loads. A
    // thumbnail that is a picture of a figure is a different claim from the
    // figure.
    var list = (example.make ? example.make(opensAs(example)) : example.files);
    var answer = K.run(example.command, list, 360);
    if (answer.ok) {
      into.innerHTML = answer.body;
      unique(into, "pv" + previews++ + "-");
    } else {
      into.textContent = "";
    }
  }

  // The previews used to be drawn as the cards were built, all two dozen of
  // them, before the dialog painted anything. That was a millisecond each when
  // every example carried a few hundred bytes; with data worth looking at it
  // is closer to seventy, and the dialog took 1.7 seconds to open on a click.
  //
  // They are drawn a few at a time now, after the dialog is up, so the panel
  // appears at once and the pictures arrive into it. A card with no picture yet
  // still carries its name, its group and its flags, which is what the reader
  // is reading while they wait.
  var waiting = [];
  var draining = 0;

  function drain() {
    if (draining) clearTimeout(draining);
    draining = setTimeout(function () {
      var until = Date.now() + 12;
      while (waiting.length && Date.now() < until) {
        var next = waiting.shift();
        preview(next[0], next[1]);
      }
      draining = 0;
      if (waiting.length) drain();
    }, 0);
  }

  function cards(filter) {
    var body = el.pickerBody;
    body.textContent = "";
    waiting.length = 0;
    var needle = (filter || "").trim().toLowerCase();
    var shown = 0;
    var groups = [];
    EXAMPLES.forEach(function (example) {
      var hay = (example.name + " " + example.group + " " + example.command).toLowerCase();
      if (needle && hay.indexOf(needle) < 0) return;
      var set = groups.filter(function (g) { return g.name === example.group; })[0];
      if (!set) { set = { name: example.group, items: [] }; groups.push(set); }
      set.items.push(example);
      shown++;
    });

    if (!shown) {
      var none = document.createElement("p");
      none.className = "pg-empty";
      none.textContent = "Nothing here answers to " + JSON.stringify(filter) + ".";
      body.appendChild(none);
      return;
    }

    // Without the program there is nothing to draw the previews with, and
    // twenty-one empty boxes read as twenty-one figures that failed. The cards
    // still carry their command and their files, so they are still worth
    // opening; the boxes are what go.
    if (!K.ready()) {
      var note = document.createElement("p");
      note.className = "pg-empty";
      note.textContent =
        "The program has not arrived, so these have no preview. The commands " +
        "and their files are still here, and they run in a terminal.";
      body.appendChild(note);
    }

    groups.forEach(function (set) {
      var section = document.createElement("section");
      section.className = "pg-set";
      var head = document.createElement("h3");
      head.textContent = set.name;
      section.appendChild(head);

      var grid = document.createElement("div");
      grid.className = "pg-grid";
      set.items.forEach(function (example) {
        var card = document.createElement("button");
        card.type = "button";
        card.className = "pg-card";

        var shot = null;
        if (K.ready()) {
          shot = document.createElement("div");
          shot.className = "pg-preview";
          card.appendChild(shot);
        }

        var title = document.createElement("span");
        title.className = "pg-card-title";
        title.textContent = example.name;
        card.appendChild(title);

        var meta = document.createElement("span");
        meta.className = "pg-card-meta";
        // The flags it uses, which is what a reader is actually shopping for.
        meta.textContent = K.words(example.command)
          .filter(function (word) { return word.indexOf("--") === 0; })
          .slice(0, 3)
          .join(" ");
        card.appendChild(meta);

        card.addEventListener("click", function () {
          closePicker();
          load(example);
        });
        grid.appendChild(card);
        // Queued rather than drawn here, and drained on a timer rather than on
        // a frame. A preview needs no layout to have happened, and a tab that
        // is not on screen never gets a frame at all, so a panel opened in a
        // background tab would come up with every card empty; a timeout still
        // fires there.
        if (shot) waiting.push([example, shot]);
      });
      section.appendChild(grid);
      body.appendChild(section);
    });
    drain();
  }

  function openPicker() {
    lastFocus = document.activeElement;
    opened = true;
    el.examples.setAttribute("aria-expanded", "true");
    el.search.value = "";
    cards("");
    // `showModal` rather than an attribute, so the page behind is inert and
    // the focus cannot leave. Escape is the element's own.
    el.picker.showModal();
    el.search.focus();
  }

  function closePicker() {
    if (el.picker.open) el.picker.close();
    // Done here rather than on the element's `close` event, which does not
    // fire in every engine: a listener added in front of a `close()` call was
    // measured never running, and the button was left saying it was expanded
    // over a panel that had gone.
    opened = false;
    el.examples.setAttribute("aria-expanded", "false");
    if (lastFocus) lastFocus.focus();
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

  // Where the splitter sits, as a number, so the keyboard and the pointer move
  // the same thing and the element can say where it is.
  var at = 38;

  function setSplit(fraction) {
    at = Math.min(75, Math.max(20, fraction));
    el.panes.style.setProperty("--pg-split", at.toFixed(1) + "%");
    el.split.setAttribute("aria-valuenow", String(Math.round(at)));
    draw();
  }

  function splitKeys(event) {
    // It announced itself as a separator with a value and answered to no key
    // at all, which is a control that says it can be adjusted and cannot.
    var stacked = el.panes.classList.contains("pg-stacked");
    var less = stacked ? "ArrowUp" : "ArrowLeft";
    var more = stacked ? "ArrowDown" : "ArrowRight";
    if (event.key === less) setSplit(at - (event.shiftKey ? 10 : 2));
    else if (event.key === more) setSplit(at + (event.shiftKey ? 10 : 2));
    else if (event.key === "Home") setSplit(20);
    else if (event.key === "End") setSplit(75);
    else return;
    event.preventDefault();
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
      setSplit(fraction * 100);
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
               "full", "examples", "picker", "search", "controls"];
    ids.forEach(function (name) { el[name] = document.getElementById("pg-" + name); });
    el.pickerBody = document.getElementById("pg-picker-body");
    el.pickerClose = document.getElementById("pg-picker-close");
    if (!el.app || !el.command) return;
    el.app.hidden = false;

    el.examples.addEventListener("click", openPicker);
    el.pickerClose.addEventListener("click", closePicker);
    el.search.addEventListener("input", function () { cards(el.search.value); });
    // Clicking the dimmed page behind the panel closes it, and so does Escape,
    // which are the two things every panel like this answers to.
    el.picker.addEventListener("mousedown", function (event) {
      if (event.target === el.picker) closePicker();
    });
    // Escape, the close button and the backdrop all end at the same place.
    // The element closes itself on Escape, so that path is taken over rather
    // than left to an event that may not arrive.
    el.picker.addEventListener("cancel", function (event) {
      event.preventDefault();
      closePicker();
    });
    el.picker.addEventListener("keydown", function (event) {
      if (event.key === "Escape") { event.preventDefault(); closePicker(); }
    });
    el.picker.addEventListener("close", closePicker);

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
      // The separator's orientation is a fact about the layout, and it was
      // frozen at the one the page happened to start in.
      el.split.setAttribute(
        "aria-orientation",
        el.panes.classList.contains("pg-stacked") ? "horizontal" : "vertical"
      );
      draw();
    });
    el.export.addEventListener("click", exportSvg);
    el.full.addEventListener("click", function () {
      if (document.fullscreenElement) document.exitFullscreen();
      else if (el.app.requestFullscreen) el.app.requestFullscreen();
    });
    dragSplit();
    el.split.addEventListener("keydown", splitKeys);

    // The figure answers to the keyboard here as it does on the home page.
    // Interactive turned on pan and zoom that only a pointer could reach.
    el.plot.addEventListener("keydown", function (event) {
      if (!el.live.checked) return;
      var step = 0.1;
      var edge = bounds();
      if (event.key === "ArrowLeft") move(K.within(K.panned(el.command.value, step), edge));
      else if (event.key === "ArrowRight") move(K.within(K.panned(el.command.value, -step), edge));
      else if (event.key === "+" || event.key === "=") move(K.within(K.zoomed(el.command.value, 0.8, 0.5), edge));
      else if (event.key === "-" || event.key === "_") move(K.within(K.zoomed(el.command.value, 1.25, 0.5), edge));
      else return;
      event.preventDefault();
    });

    var wide = null;
    window.addEventListener("resize", function () {
      clearTimeout(wide);
      wide = setTimeout(draw, 150);
    });
    // The palette toggle rewrites an attribute rather than reloading, so the
    // figure has to be told.
    // The palette toggle rewrites an attribute rather than reloading, so both
    // the figure and any previews already on screen have to be told: a preview
    // drawn light and left on a dark panel is the stale one, not the new one.
    K.onScheme(function () {
      draw();
      if (opened) cards(el.search.value);
    });

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
        el.plot.removeAttribute("role");
        el.plot.removeAttribute("aria-label");
        // A message a reader can act on rather than a file path they cannot.
        el.plot.textContent =
          "The program that draws the figure did not arrive, so nothing on " +
          "this page can run. Every command here runs the same way in a " +
          "terminal: ";
        var guide = document.createElement("a");
        guide.href = "../guide/cli/";
        guide.textContent = "the command line guide";
        el.plot.appendChild(guide);
        el.plot.appendChild(document.createTextNode("."));
      });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start);
  } else {
    start();
  }
})();
