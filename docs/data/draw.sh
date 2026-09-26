#!/bin/sh
# Draws the figures the "Start here" and "Your data" pages show, from the
# files in this folder, with exactly the commands those pages print. Each is
# drawn twice, for the light page and the dark one.
#
#   sh docs/data/draw.sh path/to/karyon
set -e
karyon="${1:-karyon}"
here="$(cd "$(dirname "$0")" && pwd)"
out="$here/../assets/start"
mkdir -p "$out"
cd "$here"
# Each is drawn on the colour of the page it sits on, light or dark, so it is
# part of the page rather than a picture laid on it.
draw() {
  name="$1"; shift
  "$karyon" "$@" --width 720 2>/dev/null | sed 's/#ffffff/#fbfaff/g' > "$out/$name.svg"
  "$karyon" "$@" --width 720 --theme dark 2>/dev/null | sed 's/#120b2b/#0d0822/g' > "$out/$name-dark.svg"
}
draw reads rpoB reads.bam genes.gff3 calls.vcf.gz
draw scan 1 gwas.assoc --threshold genome-wide
draw tree tree.nwk --traits samples.tsv --columns lineage,country
draw alignment --msa aln.fasta --with-tree tree.nwk
draw assemblies asm1_chr1 assemblies.paf
draw genome NC_000962.3 sampleA.bedgraph sampleB.bedgraph
draw locus 1:661,000-861,000 gwas.assoc --ld lead.ld --threshold genome-wide \
  genetic_map.txt --label recombination
draw heatmap NC_000962.3 --heatmap depths.tsv --relative --with-tree tree.nwk --label depth
draw pairs rpoB genes.gff3 linkage.ld
draw time --frequencies lineages.tsv --phylodynamics reproduction.tsv --threshold 1
draw selection --selection fel.csv
draw signal reads.slow5
