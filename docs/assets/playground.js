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
      group: "Phylogeny",
      command: "tangle:1-1000 --no-axis --tanglegram before.nwk --against after.nwk",
      files: [
        { name: "before.nwk", body: "((a:1,b:1):1,(c:1,d:1):1);\n" },
        { name: "after.nwk", body: "((a:1,c:1):1,(b:1,d:1):1);\n" },
      ],
    },
    {
      name: "Gene neighbourhoods",
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
      group: "Reads and molecules",
      command: "chr1:1-160 --sequence ref.fa --label reference --pileup reads.sam --label reads",
      files: [
        { name: "ref.fa", body:
            ">chr1\n" +
            "ACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAAG\n" +
            "GACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTAA\n" +
            "GGACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTTA\n" +
            "AGGACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCTT\n" +
            "AAGGACGTTGCAACGTATGCCGATTACGGCATGCATTAGCCGGATCGATCGTTAAGGCCT\n" },
        { name: "reads.sam", body: "" },
      ],
      make: function () {
        var ref = reference();
        var sam = "";
        for (var i = 0; i < 16; i++) {
          sam += "r" + i + "\t" + (i % 2 === 0 ? 0 : 16) + "\tchr1\t" + (1 + i * 5) +
                 "\t60\t60M\t*\t0\t0\t" + ref.slice(i * 5, i * 5 + 60) + "\t*\n";
        }
        return [
          { name: "ref.fa", body: fastaOf("chr1", ref) },
          { name: "reads.sam", body: sam },
        ];
      },
    },
    {
      name: "One molecule in pieces",
      group: "Reads and molecules",
      command: "chr1:1-9,000 --split-reads split.sam --label 'split reads'",
      files: [
        { name: "split.sam", body:
            "m1\t0\tchr1\t540\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6540,-,900S700M0S,60,0;\n" +
            "m2\t0\tchr1\t580\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6580,+,900S700M0S,60,0;\n" +
            "m3\t0\tchr1\t620\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6620,-,900S700M0S,60,0;\n" +
            "m4\t0\tchr1\t660\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6660,+,900S700M0S,60,0;\n" +
            "m5\t0\tchr1\t700\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6700,-,900S700M0S,60,0;\n" +
            "m6\t0\tchr1\t740\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6740,+,900S700M0S,60,0;\n" +
            "m7\t0\tchr1\t780\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6780,-,900S700M0S,60,0;\n" +
            "m8\t0\tchr1\t820\t60\t0S900M700S\t*\t0\t0\t*\t*\tSA:Z:chr1,6820,+,900S700M0S,60,0;\n" },
      ],
    },
    {
      name: "Modified bases",
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
      group: "Phylogeny",
      command: "tree:1-100 --no-axis --tree tree.nwk --label phylogeny",
      files: [
        { name: "tree.nwk", body:
            "(((s1:0.01,s2:0.012):0.02,(s3:0.008,s4:0.011):0.018):0.03,(s5:0.02,s6:0.017):0.025);\n" },
      ],
    },
    {
      name: "Recombination on a tree",
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
      group: "Alignments and rearrangements",
      command: "chr1:1-200,000 --structural sv.vcf --label 'structural variants'",
      files: [
        { name: "sv.vcf", body:
            "chr1\t10000\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=19000;SVLEN=9000\n" +
            "chr1\t38000\t.\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=50000;SVLEN=12000\n" +
            "chr1\t66000\t.\tN\t<INV>\t.\t.\tSVTYPE=INV;END=73000;SVLEN=7000\n" +
            "chr1\t94000\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=109000;SVLEN=15000\n" +
            "chr1\t122000\t.\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=128000;SVLEN=6000\n" +
            "chr1\t150000\t.\tN\t<INV>\t.\t.\tSVTYPE=INV;END=161000;SVLEN=11000\n" },
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
      // Rebuilt rather than reassigned: writing the whole attribute took
      // `pg-dragging` off one pixel into a drag, so the grabbing cursor
      // flickered back to a hand for the rest of it.
      el.plot.classList.remove("pg-failed");
      el.plot.classList.add("pg-plot");
      el.plot.classList.toggle("pg-live", el.live.checked);
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
    // Only once the figure has focus, so a wheel over a figure nobody has
    // clicked scrolls the page rather than being eaten by it.
    if (!el.live.checked || document.activeElement !== el.plot) return;
    event.preventDefault();
    var box = el.plot.getBoundingClientRect();
    var at = Math.min(1, Math.max(0, (event.clientX - box.left) / box.width));
    move(K.zoomed(el.command.value, event.deltaY > 0 ? 1.25 : 0.8, at));
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

  function preview(example, into) {
    // Drawn by the program, out of the same files the example loads. A
    // thumbnail that is a picture of a figure is a different claim from the
    // figure, and at a millisecond each there is no reason to make it.
    var list = (example.make ? example.make() : example.files);
    var answer = K.run(example.command, list, 360);
    if (answer.ok) {
      into.innerHTML = answer.body;
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
               "full", "examples", "picker", "search"];
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
      if (event.key === "ArrowLeft") move(K.panned(el.command.value, step));
      else if (event.key === "ArrowRight") move(K.panned(el.command.value, -step));
      else if (event.key === "+" || event.key === "=") move(K.zoomed(el.command.value, 0.8, 0.5));
      else if (event.key === "-" || event.key === "_") move(K.zoomed(el.command.value, 1.25, 0.5));
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
