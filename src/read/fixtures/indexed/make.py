#!/usr/bin/env python3
"""Writes the bgzipped text files the command line's index tests read, and the
indexes htslib writes for them.

Each file is made up here from a seed, so the same run writes the same rows,
and then compressed and indexed by the tools a reader's own files come from:
bgzip and tabix, and bcftools for a VCF's .csi. The tests read each file a
window at a time through its index and whole, and compare what they draw, so
they hold the reader to what tabix wrote, not to what this script thought it
should write.

The tools are htslib 1.24's (bgzip, tabix, bcftools). Run from anywhere:

    python3 src/read/fixtures/indexed/make.py

bgzip writes blocks of 64 KiB, so a file of a few hundred rows would be one
block and its index one stretch. Each file is cut into blocks of a few
kilobytes instead, by compressing pieces of it and joining them as bgzip would
have written a larger file, so a window is read out of the middle of a file of
many blocks as it is out of a large one. The VCF of a cohort is cut anywhere,
rows across two blocks included, as `bgzip --binary` and bgzip before 1.16
cut; the others where a line ends, as bgzip cuts by default now.
"""

import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# The empty block bgzip ends a file with, left off every piece but the last.
EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzip(text, name, piece=3000, anywhere=False):
    """Compresses `text` into `name` in blocks of about `piece` bytes."""
    data = text.encode()
    out = bytearray()
    at = 0
    while at < len(data):
        end = min(at + piece, len(data))
        if not anywhere and end < len(data):
            end = data.index(b"\n", end) + 1
        block = subprocess.run(
            ["bgzip", "-c"], input=data[at:end], capture_output=True, check=True
        ).stdout
        assert block.endswith(EOF)
        out += block[: -len(EOF)]
        at = end
    out += EOF
    with open(os.path.join(HERE, name), "wb") as file:
        file.write(out)


def run(*command):
    subprocess.run(command, cwd=HERE, check=True)


def cohort(rng):
    """Calls of three samples on three sequences: substitutions, insertions,
    deletions of up to 300 bases, sites of two alternates, and a fraction on
    most of them."""
    lengths = {"chr1": 3_000_000, "chr2": 1_000_000, "chr3": 200_000}
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={name},length={length}>" for name, length in lengths.items()]
    lines += [
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    ]
    for name, rows in (("chr1", 300), ("chr2", 120), ("chr3", 30)):
        for pos in sorted(rng.sample(range(1, lengths[name] - 700), rows)):
            ref = rng.choice("ACGT")
            shape = rng.random()
            if shape < 0.15:
                ref += "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 300)))
                alt = ref[0]
            elif shape < 0.3:
                alt = ref + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 8)))
            elif shape < 0.4:
                alt = ",".join(rng.sample([b for b in "ACGT" if b != ref], 2))
            else:
                alt = rng.choice([b for b in "ACGT" if b != ref])
            alleles = alt.count(",") + 1
            info = "."
            if rng.random() < 0.8:
                info = "AF=" + ",".join(f"{rng.random():.3f}" for _ in range(alleles))
            calls = [
                f"{rng.randint(0, alleles)}/{rng.randint(0, alleles)}"
                if rng.random() < 0.9
                else "./."
                for _ in range(3)
            ]
            lines.append(
                f"{name}\t{pos}\t{name}_{pos}\t{ref}\t{alt}\t50\tPASS\t{info}\tGT\t"
                + "\t".join(calls)
            )
    bgzip("\n".join(lines) + "\n", "cohort.vcf.gz", piece=4000, anywhere=True)
    run("tabix", "-f", "-p", "vcf", "cohort.vcf.gz")
    run("bcftools", "index", "-f", "cohort.vcf.gz")


def depth(rng):
    """A bedGraph of steps on two sequences, wide and narrow."""
    lines = []
    for name, length in (("chr1", 2_500_000), ("chr2", 600_000)):
        at = rng.randint(0, 5000)
        while at < length:
            width = rng.choice([20, 150, 900, 4000, 30_000])
            lines.append(f"{name}\t{at}\t{at + width}\t{rng.randint(0, 80)}")
            at += width + rng.choice([0, 0, 0, 50, 7000])
    bgzip("\n".join(lines) + "\n", "depth.bedgraph.gz")
    run("tabix", "-f", "-p", "bed", "depth.bedgraph.gz")


def genes(rng):
    """Genes of one to three transcripts and two to eight exons, as GFF3 with
    the region row NCBI opens a sequence with, and as BED12."""
    gff = ["##gff-version 3"]
    bed = []
    rows = []
    for name, length in (("chr1", 2_000_000), ("chr2", 500_000)):
        gff.append(f"##sequence-region {name} 1 {length}")
        rows.append((name, 1, f"{name}\tRefSeq\tregion\t1\t{length}\t.\t+\t.\tID={name}:1..{length}"))
        at = rng.randint(1000, 20_000)
        number = 0
        while at < length - 60_000:
            number += 1
            gene = f"{name}g{number}"
            strand = rng.choice("+-")
            exons = []
            cursor = at
            for _ in range(rng.randint(2, 8)):
                size = rng.randint(80, 1500)
                exons.append((cursor, cursor + size - 1))
                cursor += size + rng.randint(200, 9000)
            end = exons[-1][1]
            rows.append((name, at, f"{name}\t.\tgene\t{at}\t{end}\t.\t{strand}\t.\tID={gene};Name={gene}"))
            for t in range(rng.randint(1, 3)):
                kept = [exons[0]] + [e for e in exons[1:-1] if rng.random() < 0.7] + [exons[-1]]
                transcript = f"{gene}.t{t + 1}"
                rows.append(
                    (name, kept[0][0], f"{name}\t.\tmRNA\t{kept[0][0]}\t{kept[-1][1]}\t.\t{strand}\t.\tID={transcript};Parent={gene}")
                )
                for start, stop in kept:
                    rows.append((name, start, f"{name}\t.\texon\t{start}\t{stop}\t.\t{strand}\t.\tParent={transcript}"))
                    rows.append((name, start, f"{name}\t.\tCDS\t{start}\t{stop}\t.\t{strand}\t0\tParent={transcript}"))
            sizes = ",".join(str(stop - start + 1) for start, stop in exons)
            starts = ",".join(str(start - at) for start, _ in exons)
            bed.append(
                f"{name}\t{at - 1}\t{end}\t{gene}\t0\t{strand}\t{at - 1}\t{end}\t0\t{len(exons)}\t{sizes},\t{starts},"
            )
            at = end + rng.randint(-20_000, 60_000)
            at = max(at, exons[0][0] + 1)
    # Sorted as `sort -k1,1 -k4,4n` sorts it, which tabix asks for.
    rows.sort(key=lambda row: (row[0], row[1]))
    gff += [row[2] for row in rows]
    bgzip("\n".join(gff) + "\n", "genes.gff3.gz")
    run("tabix", "-f", "-p", "gff", "genes.gff3.gz")
    run("tabix", "-f", "-C", "-p", "gff", "genes.gff3.gz")
    bed.sort(key=lambda row: (row.split("\t")[0], int(row.split("\t")[1])))
    bgzip("\n".join(bed) + "\n", "genes.bed.gz")
    run("tabix", "-f", "-p", "bed", "genes.bed.gz")
    # The same genes as UCSC's GTF writes them, exons alone.
    gtf = []
    for row in rows:
        cols = row[2].split("\t")
        if cols[2] != "exon":
            continue
        transcript = cols[8].split("=")[1]
        gene = transcript.split(".")[0]
        cols[8] = f'gene_id "{gene}"; transcript_id "{transcript}";'
        gtf.append("\t".join(cols))
    bgzip("\n".join(gtf) + "\n", "exons.gtf.gz")
    run("tabix", "-f", "-p", "gff", "exons.gtf.gz")


def methylation(rng):
    """modkit's bedMethyl of a dual-mode run: m and h at every cytosine, on
    both strands."""
    lines = []
    for name, length in (("chr1", 1_500_000), ("chr2", 300_000)):
        for start in sorted(rng.sample(range(0, length), 150)):
            strand = rng.choice("+-")
            for code in "hm":
                valid = rng.randint(0, 40)
                modified = rng.randint(0, valid)
                fraction = f"{100 * modified / valid:.2f}" if valid else "0.00"
                lines.append(
                    f"{name}\t{start}\t{start + 1}\t{code}\t{valid}\t{strand}\t{start}\t{start + 1}\t"
                    f"255,0,0\t{valid}\t{fraction}\t{modified}\t{valid - modified}\t0\t0\t0\t0\t0"
                )
    bgzip("\n".join(lines) + "\n", "methyl.bed.gz")
    run("tabix", "-f", "-p", "bed", "methyl.bed.gz")


def junctions(rng):
    """STAR's SJ.out.tab: introns counted from one, inclusive."""
    lines = []
    for name, length in (("chr1", 1_500_000), ("chr2", 400_000)):
        for first in sorted(rng.sample(range(1, length - 20_000), 250)):
            last = first + rng.randint(60, 15_000)
            lines.append(
                f"{name}\t{first}\t{last}\t{rng.randint(0, 2)}\t{rng.randint(0, 6)}\t{rng.randint(0, 1)}\t"
                f"{rng.choice([0, rng.randint(1, 90)])}\t{rng.randint(0, 9)}\t{rng.randint(5, 60)}"
            )
    bgzip("\n".join(lines) + "\n", "sj.tab.gz")
    run("tabix", "-f", "-s1", "-b2", "-e3", "sj.tab.gz")


def scans(rng):
    """An association scan as PLINK 2 writes it, its header behind a `#`, and
    as PLINK 1 does, with a plain header, which tabix is told to skip."""
    two = ["#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP"]
    one = ["CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP"]
    for name, length in (("1", 2_000_000), ("2", 800_000)):
        for pos in sorted(rng.sample(range(1, length), 250)):
            p = 10 ** -rng.uniform(0, 9)
            beta = rng.uniform(-1, 1)
            two.append(f"{name}\t{pos}\trs{pos}\tA\tG\tG\tADD\t500\t{beta:.4f}\t0.1\t{beta * 10:.3f}\t{p:.4g}")
            one.append(f"{name}\trs{pos}\t{pos}\tG\tADD\t500\t{beta:.4f}\t{beta * 10:.3f}\t{p:.4g}")
    bgzip("\n".join(two) + "\n", "scan2.tsv.gz")
    run("tabix", "-f", "-s1", "-b2", "-e2", "scan2.tsv.gz")
    bgzip("\n".join(one) + "\n", "scan1.tsv.gz")
    run("tabix", "-f", "-S1", "-s1", "-b3", "-e3", "scan1.tsv.gz")


def windows(rng):
    """A table of windows, a value per sample, as `bedtools unionbedg -header`
    writes it, 0-based and with a plain header."""
    lines = ["chrom\tstart\tend\tS1\tS2\tS3"]
    for name, length in (("chr1", 2_000_000), ("chr2", 500_000)):
        for start in range(0, length, 5000):
            values = "\t".join(f"{rng.uniform(0, 60):.2f}" if rng.random() < 0.95 else "NA" for _ in range(3))
            lines.append(f"{name}\t{start}\t{start + 5000}\t{values}")
    bgzip("\n".join(lines) + "\n", "windows.tsv.gz")
    run("tabix", "-f", "-0", "-S1", "-s1", "-b2", "-e3", "windows.tsv.gz")


def read_whole(rng):
    """Three files with an index that are read whole all the same, since the
    rows over a window do not say what the whole file says of them: a scan
    with no header, whose values are p-values only if every one of them lies
    between 0 and 1; a table of windows in the long form, a sample to a row,
    which names its samples on its rows and leaves one out of some windows;
    and a BED whose first row has a dot in column seven, where GFF3 writes a
    strand, so the whole file is taken for GFF3 and its rows over a window
    further on for BED."""
    lines = []
    for name, length in (("1", 1_000_000), ("2", 300_000)):
        for pos in sorted(rng.sample(range(1, length), 120)):
            # Above 1 only far along the first sequence.
            value = rng.uniform(1, 9) if name == "1" and pos > 800_000 else rng.random()
            lines.append(f"{name}\t{pos}\t{value:.4g}")
    bgzip("\n".join(lines) + "\n", "scan0.tsv.gz", piece=1000)
    run("tabix", "-f", "-s1", "-b2", "-e2", "scan0.tsv.gz")
    lines = ["chrom\tstart\tend\tsample\tdepth"]
    for start in range(0, 400_000, 10_000):
        for sample in ("S1", "S2", "S3"):
            if sample == "S3" and start < 200_000:
                continue
            lines.append(f"chr1\t{start}\t{start + 10_000}\t{sample}\t{rng.uniform(0, 50):.2f}")
    bgzip("\n".join(lines) + "\n", "long.tsv.gz", piece=1000)
    run("tabix", "-f", "-0", "-S1", "-s1", "-b2", "-e3", "long.tsv.gz")
    lines = ["chr1\t100\t900\tfirst\t0\t+\t.\t."]
    for start in range(10_000, 400_000, 2_500):
        lines.append(f"chr1\t{start}\t{start + 800}\tg{start}\t0\t+\t{start}\t{start + 800}")
    bgzip("\n".join(lines) + "\n", "odd.bed.gz", piece=1000)
    run("tabix", "-f", "-p", "bed", "odd.bed.gz")


def structural():
    """Breakends joined across 800 kb, which the lower row draws as one arc,
    and a deletion with an END."""
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t3000\tbnd1\tN\tN[chr1:800000[\t50\tPASS\tSVTYPE=BND;MATEID=bnd2",
        "chr1\t400500\tdel1\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=401500",
        "chr1\t800000\tbnd2\tN\t]chr1:3000]N\t50\tPASS\tSVTYPE=BND;MATEID=bnd1",
    ]
    bgzip("\n".join(lines) + "\n", "sv.vcf.gz", piece=80)
    run("tabix", "-f", "-p", "vcf", "sv.vcf.gz")


def seven():
    """Calls whose row on chr2 has seven columns, which tabix indexes, since
    it reads no further than REF, and the reader of calls refuses."""
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t100\tv1\tA\tG\t50\tPASS\t.",
        "chr1\t250\tv2\tC\tT\t50\tPASS\t.",
        "chr2\t300\tv3\tG\tA\t50\tPASS",
    ]
    bgzip("\n".join(lines) + "\n", "seven.vcf.gz", piece=60)
    run("tabix", "-f", "-p", "vcf", "seven.vcf.gz")


def main():
    for tool in ("bgzip", "tabix", "bcftools"):
        version = subprocess.run([tool, "--version"], capture_output=True, text=True).stdout
        print(version.splitlines()[0], file=sys.stderr)
    cohort(random.Random(1))
    depth(random.Random(2))
    genes(random.Random(3))
    methylation(random.Random(4))
    junctions(random.Random(5))
    scans(random.Random(6))
    windows(random.Random(7))
    read_whole(random.Random(8))
    structural()
    seven()


if __name__ == "__main__":
    main()
