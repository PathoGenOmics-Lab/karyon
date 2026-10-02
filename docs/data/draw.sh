#!/bin/sh
# Draws the figures the "Start here" and "Your data" pages show, from the
# files in this folder, with exactly the commands those pages print. Each is
# drawn twice, for the light page and the dark one.
#
#   sh docs/data/draw.sh path/to/karyon
#   sh docs/data/draw.sh path/to/karyon pdf some/folder
#
# The second form draws the same figures with the same commands as PDF, into
# a folder of its own, for CI to read back with Ghostscript and poppler. It
# goes through the command line's own `-o`, which is the way a reader of the
# site gets one; no PDF is committed.
set -e
karyon="${1:-karyon}"
kind="${2:-svg}"
here="$(cd "$(dirname "$0")" && pwd)"
out="$here/../assets/start"
if [ "$kind" = pdf ]; then
  mkdir -p "${3:?a folder to draw the PDFs into}"
  out="$(cd "$3" && pwd)"
fi
mkdir -p "$out"
cd "$here"
# Each is drawn on the colour of the page it sits on, light or dark, so it is
# part of the page rather than a picture laid on it, and is what the program
# in the page draws from the same command once it arrives: a test in the
# playground holds the two to one figure.
draw() {
  name="$1"; shift
  if [ "$kind" = pdf ]; then
    "$karyon" "$@" --width 720 --background '#fbfaff' -o "$out/$name.pdf"
    "$karyon" "$@" --width 720 --theme dark --background '#0d0822' -o "$out/$name-dark.pdf"
    return
  fi
  "$karyon" "$@" --width 720 --background '#fbfaff' 2>/dev/null > "$out/$name.svg"
  "$karyon" "$@" --width 720 --theme dark --background '#0d0822' 2>/dev/null > "$out/$name-dark.svg"
}
draw reads rpoB reads.bam genes.gff3 calls.vcf.gz
draw scan 1 gwas.assoc --threshold genome-wide
draw genome-scan trait.assoc --threshold genome-wide
draw tree tree.nwk --traits samples.tsv --columns lineage,country
draw alignment --msa aln.fasta --with-tree tree.nwk
draw assemblies asm1_chr1 assemblies.paf
draw genome NC_000962.3 sampleA.bedgraph sampleB.bedgraph --same-scale
draw locus 1:661,000-861,000 gwas.assoc --ld lead.ld --threshold genome-wide \
  --with-recombination genetic_map.txt
draw heatmap NC_000962.3 --heatmap depths.tsv --relative --with-tree tree.nwk --label depth
draw pairs rpoB genes.gff3 linkage.ld
draw time --frequencies lineages.tsv --phylodynamics reproduction.tsv --threshold 1
draw selection --selection fel.csv
draw signal reads.slow5
