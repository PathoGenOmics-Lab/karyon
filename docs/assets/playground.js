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
  function reference() {
    var unit = "ACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAGG";
    var out = "";
    for (var i = 0; i < 300; i++) out += unit[i % 61];
    return out;
  }

  function fastaOf(name, seq) {
    var out = ">" + name + "\n";
    for (var i = 0; i < seq.length; i += 60) out += seq.slice(i, i + 60) + "\n";
    return out;
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
        for (var at = 757000; at < 767000; at += 100) {
          var dip = at > 761890 && at < 762030 ? 0.06 : 1;
          depth += "NC_000962.3 " + at + " " + (at + 100) + " " +
                   Math.round((56 + 8 * Math.sin((at - 757000) / 900)) * dip) + "\n";
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
        { kind: "choice", flag: "--aggregate", after: "--coverage", label: "How a pixel that covers many bases chooses", options: ["max", "mean", "min"] },
        { kind: "toggle", flag: "--log", after: "--coverage", label: "A log scale for the depth" },
      ],
      group: "Signal and annotation",
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
      bounds: { from: 1, to: 1000, min: 60 },
      controls: [
        { kind: "note", label: "Two phylogenies face to face. The axis is not a coordinate, so panning it would mean nothing; what changes a tanglegram is which trees you give it." },
        { kind: "toggle", flag: "--no-axis", label: "The coordinate ruler, which measures nothing here" },
      ],
      group: "Phylogeny",
      command: "tangle:1-1000 --no-axis --tanglegram before.nwk --against after.nwk",
      files: [
        { name: "before.nwk", body: "((a:1,b:1):1,(c:1,d:1):1);\n" },
        { name: "after.nwk", body: "((a:1,c:1):1,(b:1,d:1):1);\n" },
      ],
    },
    {
      name: "Gene neighbourhoods",
      bounds: { from: 1, to: 4000, min: 1001 },
      controls: [
        { kind: "region" },
      ],
      group: "Comparisons across genomes",
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
      bounds: { from: 1, to: 700, min: 233 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--analysis", after: "--domains", label: "Which member database annotated it", options: ["Pfam"] },
      ],
      group: "Comparisons across genomes",
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
      bounds: { from: 1, to: 200, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--context", after: "--bisulfite", label: "Which cytosine context", options: ["CpG"] },
      ],
      group: "Reads and molecules",
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
    {
      name: "Reference bases",
      bounds: { from: 1, to: 300, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Signal and annotation",
      command: "chr1:1-120 --axis --sequence ref.fa --label reference --orfs ref.fa --label 'reading frames'",
      files: [
        { name: "ref.fa", body:
            ">chr1\n" +
            "ACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "GACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAA\n" +
            "GGACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTA\n" +
            "AGGACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTT\n" +
            "AAGGACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCT\n" },
      ],
    },
    {
      name: "A cytogenetic ideogram",
      bounds: { from: 1, to: 2000000, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Signal and annotation",
      command: "chr1:1-2,000,000 --ideogram bands.txt --label chromosome",
      files: [
        { name: "bands.txt", body:
            "chr1 0 200000 p11 gneg\n" +
            "chr1 200000 400000 p10 gpos25\n" +
            "chr1 400000 600000 p9 gpos50\n" +
            "chr1 600000 800000 p8 gpos75\n" +
            "chr1 800000 1000000 p7 gpos100\n" +
            "chr1 1000000 1200000 q1 acen\n" +
            "chr1 1200000 1400000 q2 gvar\n" +
            "chr1 1400000 1600000 q3 stalk\n" +
            "chr1 1600000 1800000 q4 gneg\n" +
            "chr1 1800000 2000000 q5 gpos50\n" },
      ],
    },
    {
      name: "A genome-wide scan",
      bounds: { from: 1, to: 1000000, min: 1997 },
      controls: [
        { kind: "region" },
      ],
      group: "Association and genotype",
      command: "chr1:1-1,000,000 --manhattan assoc.tsv --label association",
      files: [
        { name: "assoc.tsv", body: "" },
      ],
      // Five hundred rows built here rather than pasted, so the example is one
      // a reader can pan across rather than one they can read.
      make: function () {
        var rows = "pos\tp\n";
        for (var i = 0; i < 500; i++) {
          var peak = 6 * Math.exp(-Math.pow(i - 250, 2) / 900);
          var p = Math.pow(10, -(1 + peak + 0.6 * Math.sin(i)));
          rows += (i * 2000 + 1) + "\t" + p.toPrecision(4) + "\n";
        }
        return [{ name: "assoc.tsv", body: rows }];
      },
    },
    {
      name: "A genotype matrix",
      bounds: { from: 1, to: 400, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Association and genotype",
      command: "chr1:1-400 --matrix geno.tsv --label 'allele fraction'",
      files: [
        { name: "geno.tsv", body:
            "sample\t51\t97\t149\t203\t258\t311\t355\n" +
            "S01\t0.95\t0.05\t0.95\t0.05\t0.95\t0.05\t0.95\n" +
            "S02\t0.95\t0.05\t0.95\t0.05\t0.95\t0.05\t0.95\n" +
            "S03\t0.05\t0.95\t0.05\t0.95\t0.05\t0.95\t0.05\n" +
            "S04\t0.05\t0.95\t0.05\t0.95\t0.05\t0.95\t0.05\n" +
            "S05\t0.05\t0.95\t0.05\t0.95\t0.05\t0.95\tNA\n" +
            "S06\t0.95\t0.05\t0.95\t0.05\t0.95\tNA\t0.95\n" +
            "S07\t0.95\t0.05\t0.95\t0.05\tNA\t0.05\t0.95\n" +
            "S08\t0.95\t0.05\t0.95\tNA\t0.95\t0.05\t0.95\n" +
            "S09\t0.05\t0.95\tNA\t0.95\t0.05\t0.95\t0.05\n" +
            "S10\t0.05\tNA\t0.05\t0.95\t0.05\t0.95\t0.05\n" +
            "S11\tNA\t0.95\t0.05\t0.95\t0.05\t0.95\t0.05\n" +
            "S12\t0.95\t0.05\t0.95\t0.05\t0.95\t0.05\t0.95\n" },
      ],
    },
    {
      name: "A read pileup",
      bounds: { from: 1, to: 300, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Reads and molecules",
      command: "chr1:1-300 --sequence ref.fa --label reference --pileup reads.sam --label reads",
      files: [
        { name: "ref.fa", body: "" },
        { name: "reads.sam", body: "" },
      ],
      make: function () {
        var ref = reference();
        var sam = "";
        var n = 0;
        for (var at = 1; at < 260; at += 6) {
          for (var copy = 0; copy < 3; copy++) {
            var len = 55 + ((at + copy * 7) % 12);
            sam += "r" + (n++) + "\t" + (copy % 2 === 0 ? 0 : 16) + "\tchr1\t" + at +
                   "\t60\t" + len + "M\t*\t0\t0\t" + ref.slice(at - 1, at - 1 + len) + "\t*\n";
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
      bounds: { from: 1, to: 9000, min: 200 },
      controls: [
        { kind: "region" },
      ],
      group: "Reads and molecules",
      command: "chr1:1-9,000 --split-reads split.sam --label 'split reads'",
      files: [{ name: "split.sam", body: "" }],
      make: function () {
        var sam = "";
        for (var i = 0; i < 30; i++) {
          var a = 200 + i * 280;
          var b = a + 900 + (i % 5) * 120;
          var back = i % 3 === 0;
          sam += "m" + i + "\t0\tchr1\t" + a + "\t60\t600M700S\t*\t0\t0\t*\t*\t" +
                 "SA:Z:chr1," + b + "," + (back ? "-" : "+") + ",600S700M,60,0;\n";
        }
        return [{ name: "split.sam", body: sam }];
      },
    },
    {
      name: "Modified bases",
      bounds: { from: 1, to: 1000, min: 60 },
      controls: [
        { kind: "region" },
        { kind: "choice", flag: "--modification", after: "--methylation", label: "Which modified base to draw", options: ["m"] },
      ],
      group: "Reads and molecules",
      command: "chr1:1-1,000 --methylation calls.bed --modification m --label '5mC'",
      files: [
        { name: "calls.bed", body: "" },
      ],
      make: function () {
        var rows = "";
        for (var p = 20; p < 980; p += 24) {
          var pct = p < 400 ? 90 : 12;
          var mod = Math.round(40 * pct / 100);
          rows += "chr1 " + p + " " + (p + 1) + " m 40 + " + p + " " + (p + 1) +
                  " 0,0,0 40 " + pct.toFixed(2) + " " + mod + " " + (40 - mod) +
                  " 0 0 0 0 0\n";
        }
        return [{ name: "calls.bed", body: rows }];
      },
    },
    {
      name: "Alignment ribbons",
      bounds: { from: 1, to: 40000, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Alignments and rearrangements",
      command: "ctg1:1-40,000 --synteny aln.paf --label 'against chrA'",
      files: [
        { name: "aln.paf", body:
            "ctg1\t40000\t0\t3600\t-\tchrA\t60000\t5000\t8600\t3400\t3600\t60\n" +
            "ctg1\t40000\t4000\t7600\t+\tchrA\t60000\t9200\t12800\t3400\t3600\t60\n" +
            "ctg1\t40000\t8000\t11600\t+\tchrA\t60000\t13400\t17000\t3400\t3600\t60\n" +
            "ctg1\t40000\t12000\t15600\t-\tchrA\t60000\t17600\t21200\t3400\t3600\t60\n" +
            "ctg1\t40000\t16000\t19600\t+\tchrA\t60000\t21800\t25400\t3400\t3600\t60\n" +
            "ctg1\t40000\t20000\t23600\t+\tchrA\t60000\t26000\t29600\t3400\t3600\t60\n" +
            "ctg1\t40000\t24000\t27600\t-\tchrA\t60000\t30200\t33800\t3400\t3600\t60\n" +
            "ctg1\t40000\t28000\t31600\t+\tchrA\t60000\t34400\t38000\t3400\t3600\t60\n" +
            "ctg1\t40000\t32000\t35600\t+\tchrA\t60000\t38600\t42200\t3400\t3600\t60\n" },
      ],
    },
    {
      name: "The same PAF as a dot plot",
      bounds: { from: 1, to: 40000, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Alignments and rearrangements",
      command: "ctg1:1-40,000 --dotplot aln.paf --label 'dot plot'",
      files: [
        { name: "aln.paf", body:
            "ctg1\t40000\t0\t3600\t-\tchrA\t60000\t5000\t8600\t3400\t3600\t60\n" +
            "ctg1\t40000\t4000\t7600\t+\tchrA\t60000\t9200\t12800\t3400\t3600\t60\n" +
            "ctg1\t40000\t8000\t11600\t+\tchrA\t60000\t13400\t17000\t3400\t3600\t60\n" +
            "ctg1\t40000\t12000\t15600\t-\tchrA\t60000\t17600\t21200\t3400\t3600\t60\n" +
            "ctg1\t40000\t16000\t19600\t+\tchrA\t60000\t21800\t25400\t3400\t3600\t60\n" +
            "ctg1\t40000\t20000\t23600\t+\tchrA\t60000\t26000\t29600\t3400\t3600\t60\n" +
            "ctg1\t40000\t24000\t27600\t-\tchrA\t60000\t30200\t33800\t3400\t3600\t60\n" +
            "ctg1\t40000\t28000\t31600\t+\tchrA\t60000\t34400\t38000\t3400\t3600\t60\n" +
            "ctg1\t40000\t32000\t35600\t+\tchrA\t60000\t38600\t42200\t3400\t3600\t60\n" },
      ],
    },
    {
      name: "A multiple alignment",
      bounds: { from: 1, to: 61, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Sequence alignment",
      command: "aln:1-61 --msa aln.fa --label alignment",
      files: [
        { name: "aln.fa", body:
            ">sample1\n" +
            "ACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample2\n" +
            "ACGTTGCAACTTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample3\n" +
            "ACGTTGCAACTTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample4\n" +
            "ACGTTGCAACGTATGCCGATGACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" },
      ],
    },
    {
      name: "Only the variable sites",
      bounds: { from: 1, to: 61, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Sequence alignment",
      command: "aln:1-61 --snps aln.fa --label 'variable sites'",
      files: [
        { name: "aln.fa", body:
            ">sample1\n" +
            "ACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample2\n" +
            "ACGTTGCAACTTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample3\n" +
            "ACGTTGCAACTTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample4\n" +
            "ACGTTGCAACGTATGCCGATGACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" },
      ],
    },
    {
      name: "A sequence logo",
      bounds: { from: 1, to: 61, min: 60 },
      controls: [
        { kind: "region" },
      ],
      group: "Sequence alignment",
      command: "aln:1-61 --logo aln.fa --label logo",
      files: [
        { name: "aln.fa", body:
            ">sample1\n" +
            "ACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample2\n" +
            "ACGTTGCAACTTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample3\n" +
            "ACGTTGCAACTTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" +
            ">sample4\n" +
            "ACGTTGCAACGTATGCCGATGACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "G\n" },
      ],
    },
    {
      name: "A phylogeny",
      bounds: { from: 1, to: 100, min: 60 },
      controls: [
        { kind: "note", label: "A phylogeny has no coordinates, so there is nothing to pan across. What it has is a shape." },
        { kind: "toggle", flag: "--no-axis", label: "The coordinate ruler, which measures nothing here" },
      ],
      group: "Phylogeny",
      command: "tree:1-100 --no-axis --tree tree.nwk --label phylogeny",
      files: [
        { name: "tree.nwk", body:
            "(((s1:0.01,s2:0.012):0.02,(s3:0.008,s4:0.011):0.018):0.03,(s5:0.02,s6:0.017):0.025);\n" },
      ],
    },
    {
      name: "Recombination on a tree",
      bounds: { from: 1, to: 8000, min: 1701 },
      controls: [
        { kind: "region" },
      ],
      group: "Phylogeny",
      command: "NC_011900.1:1-8,000 --clades gubbins.gff --with-tree tree.nwk --label recombination",
      files: [
        { name: "gubbins.gff", body:
            "##gff-version 3\n" +
            "SEQUENCE\tGUBBINS\tCDS\t500\t1500\t0.000\t.\t0\tnode=\"N7\";taxa=\"s1 s2 s3\";\n" +
            "SEQUENCE\tGUBBINS\tCDS\t3000\t4200\t0.000\t.\t0\tnode=\"N2\";taxa=\"s5 s6\";\n" +
            "SEQUENCE\tGUBBINS\tCDS\t6000\t6800\t0.000\t.\t0\tnode=\"N9\";taxa=\"s3 s4\";\n" },
        { name: "tree.nwk", body:
            "(((s1:0.01,s2:0.012):0.02,(s3:0.008,s4:0.011):0.018):0.03,(s5:0.02,s6:0.017):0.025);\n" },
      ],
    },
    {
      name: "Structural variants",
      bounds: { from: 1, to: 200000, min: 2097 },
      controls: [
        { kind: "region" },
      ],
      group: "Alignments and rearrangements",
      command: "chr1:1-200,000 --structural sv.vcf --label 'structural variants'",
      files: [{ name: "sv.vcf", body: "" }],
      make: function () {
        var kinds = ["DEL", "DUP", "INV"];
        var vcf = "";
        for (var i = 0; i < 40; i++) {
          var at = 2000 + i * 4900;
          var len = 1200 + (i % 7) * 900;
          var k = kinds[i % 3];
          vcf += "chr1\t" + at + "\t.\tN\t<" + k + ">\t.\t.\tSVTYPE=" + k +
                 ";END=" + (at + len) + ";SVLEN=" + len + "\n";
        }
        return [{ name: "sv.vcf", body: vcf }];
      },
    },
    {
      name: "Segmented copy number",
      bounds: { from: 1, to: 2000000, min: 2097 },
      controls: [
        { kind: "region" },
        { kind: "note", label: "The ploidy is not in the file, which is why --ploidy is required rather than defaulted: a log ratio only becomes copies once you say what two copies means here." },
      ],
      group: "Signal and annotation",
      command:
        "chr8:1-2,000,000 --copy-number segments.cns --ploidy 2 \\\n" +
        "  --label 'copy number'",
      files: [
        { name: "segments.cns", body:
          "chromosome\tstart\tend\tgene\tlog2\tcn\tcn1\tcn2\n" +
          "chr8\t0\t400000\t-\t0.02\t2\t1\t1\n" +
          "chr8\t400000\t700000\tMYC\t1.70\t6\t4\t2\n" +
          "chr8\t700000\t1000000\t-\t-1.00\t1\t1\t0\n" +
          "chr8\t1000000\t1400000\t-\t0.00\t2\t1\t1\n" +
          "chr8\t1400000\t1700000\t-\tNA\tNA\tNA\tNA\n" +
          "chr8\t1700000\t2000000\t-\t0.58\t3\t2\t1\n" },
      ],
    },
    {
      name: "Per-base model attribution",
      bounds: { from: 1, to: 120, min: 40 },
      controls: [
        { kind: "region" },
        { kind: "note", label: "The bases themselves are the bars. A letter above the line is one the model leaned on, and one below is one it pulled away from, so the height is a signed score rather than a count." },
      ],
      group: "Sequence alignment",
      command:
        "promoter:1-120 --dynseq attribution.bg --with-sequence promoter.fa \\\n" +
        "  --label attribution",
      files: [
        { name: "attribution.bg", body: "" },
        { name: "promoter.fa", body: "" },
      ],
      make: function () {
        var bg = "";
        for (var i = 0; i < 120; i++) {
          var v;
          if (i >= 40 && i < 48) v = 0.9 - 0.05 * (i - 40);
          else if (i >= 70 && i < 76) v = -0.6 + 0.08 * (i - 70);
          else v = 0.06 * Math.sin(i / 5);
          bg += "promoter\t" + i + "\t" + (i + 1) + "\t" + v.toFixed(3) + "\n";
        }
        // A fixed sequence rather than a random one, so the figure is the same
        // every time the example is opened.
        var bases = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTATAAAAGCAGTGTGAATCG" +
                    "TTGCACCGATTAGGCATCAGTACCGGATTACAGCTTAAGCCGGATTCAGTACCGATTAGC";
        var fa = ">promoter\n" + bases.slice(0, 60) + "\n" + bases.slice(60) + "\n";
        return [
          { name: "attribution.bg", body: bg },
          { name: "promoter.fa", body: fa },
        ];
      },
    },
    {
      name: "Splice junctions",
      bounds: { from: 1, to: 7000, min: 200 },
      controls: [
        { kind: "region" },
        { kind: "note", label: "An arc per intron, thicker for the junctions more reads crossed. Multi-mapping reads are counted separately and never added in: a read that mapped in four places is one read." },
      ],
      group: "Reads and molecules",
      command: "chr1:1-7,000 --junctions SJ.out.tab --label junctions",
      files: [
        { name: "SJ.out.tab", body:
          "chr1\t1200\t2400\t1\t2\t1\t46\t3\t38\n" +
          "chr1\t1200\t3600\t1\t2\t0\t9\t1\t31\n" +
          "chr1\t2700\t3600\t1\t2\t1\t52\t4\t40\n" +
          "chr1\t4000\t5200\t2\t1\t1\t18\t0\t35\n" +
          "chr1\t5600\t6400\t1\t2\t1\t7\t2\t29\n" },
      ],
    },
  ];

  // ---------------------------------------------------------------------
  // The controls, which belong to the example
  // ---------------------------------------------------------------------
  //
  // Not the same for every one, because the figures are not the same thing.
  // A window to slide is what a signal over a chromosome has and a tanglegram
  // has not, and offering the tanglegram one anyway would be a control that
  // does nothing and says nothing about why. Every control is a flag, so
  // turning one rewrites a word of the command and the command stays the thing
  // that decides.

  var current = null;

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

  function load(example) {
    current = example;
    el.command.value = example.command;
    // The editor is holding the last example's file, not this one's, so it
    // has nothing to save.
    showing = null;
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

  function preview(example, into) {
    // Drawn by the program, out of the same files the example loads. A
    // thumbnail that is a picture of a figure is a different claim from the
    // figure, and at a millisecond each there is no reason to make it.
    var list = (example.make ? example.make() : example.files);
    var answer = K.run(example.command, list, 360);
    if (answer.ok) {
      into.innerHTML = answer.body;
      unique(into, "pv" + previews++ + "-");
    } else {
      into.textContent = "";
    }
  }

  function cards(filter) {
    var body = el.pickerBody;
    body.textContent = "";
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
        // Drawn now rather than on the next frame. A preview is drawn at a
        // fixed width and needs no layout to have happened, and a tab that is
        // not on screen never gets a frame at all, so a panel opened in a
        // background tab would have come up with every card empty.
        if (shot) preview(example, shot);
      });
      section.appendChild(grid);
      body.appendChild(section);
    });
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
