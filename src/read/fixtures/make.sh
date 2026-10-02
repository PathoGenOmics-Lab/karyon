#!/bin/sh
# Writes the bigWig, bigBed and 2bit files the readers beside this folder are
# tested on, from the text files here, with the UCSC tools that write them,
# and then reads each one back with the tool that undoes it. A reader's test
# compares what it reads with that second file, so it is held to what the
# tools themselves say the file holds, not to what this folder's author
# thought it should hold.
#
# The tools are kent 482 from bioconda (ucsc-bedgraphtobigwig and the five
# like it). Run from anywhere:
#
#     sh src/read/fixtures/make.sh
#
# The small -blockSize and -itemsPerSlot are on purpose: kent pads every node
# of an index to blockSize entries, 256 by default, so a few spans written
# that way are 13 KB of mostly zeros, and at two a node the indexes of these
# few rows already have inner nodes as well as leaves, which is what a large
# file has and a reader has to walk.
set -eu
cd "$(dirname "$0")"

bedGraphToBigWig -blockSize=2 -itemsPerSlot=2 signal.bedgraph signal.sizes signal.bw
bedGraphToBigWig -blockSize=2 -itemsPerSlot=2 -unc signal.bedgraph signal.sizes signal.unc.bw
bigWigToBedGraph signal.bw signal.bw.bedgraph

bedGraphToBigWig -blockSize=4 -itemsPerSlot=32 steps.bedgraph steps.sizes steps.bw
bigWigToBedGraph steps.bw steps.bw.bedgraph

bedToBigBed -type=bed12 -blockSize=2 -itemsPerSlot=2 genes.bed genes.sizes genes.bb
bigBedToBed genes.bb genes.bb.bed

bedToBigBed -type=bed6+4 -as=peaks.as -blockSize=2 -itemsPerSlot=2 peaks.bed genes.sizes peaks.bb
bigBedToBed peaks.bb peaks.bb.bed

# chr1's rows alone, against lengths that name every sequence: kent indexes
# only the sequences that hold data, so each of these names chr1 and no other.
awk '$1 == "chr1"' signal.bedgraph > sparse.bedgraph
bedGraphToBigWig -blockSize=2 -itemsPerSlot=2 sparse.bedgraph signal.sizes sparse.bw
awk '$1 == "chr1"' genes.bed > sparse.bed
bedToBigBed -type=bed12 -blockSize=2 -itemsPerSlot=2 sparse.bed genes.sizes sparse.bb

# Four columns whose fourth is a number, which as text is read as a signal.
bedToBigBed -type=bed4 -blockSize=2 -itemsPerSlot=2 scores.bed genes.sizes scores.bb
bigBedToBed scores.bb scores.bb.bed

faToTwoBit ref.fa ref.2bit
faToTwoBit -long ref.fa ref.long.2bit
twoBitToFa ref.2bit ref.2bit.fa
twoBitToFa -seq=chr1 -start=2 -end=50 ref.2bit ref.chr1.fa

faToTwoBit one.fa one.2bit
