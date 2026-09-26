"""Writes the example files the site's pages draw.

Every file here is synthetic. The genome coordinates, the gene names and their
positions are those of the reference the rest of the site uses, so a figure of
them looks like the real locus, but the sequence, the reads, the calls, the
scan, the trees and the samples are made up from a fixed seed. Running this
again writes the same files byte for byte.

It needs samtools, bgzip and tabix on the PATH, for the BAM and the VCF.
"""

import gzip
import os
import random
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
SEQ = "NC_000962.3"
LENGTH = 4411532
START, END = 759001, 764000  # 1-based, inclusive: the stretch the reads cover
GENES = [
    ("rpoB", 759807, 763325, "+"),
    ("rpoC", 763370, 767320, "+"),
]
# Calls inside rpoB, with the fraction of reads that carry them.
CALLS = [
    (760314, "C", "T", 1.00, "synonymous_variant"),
    (761110, "A", "T", 0.35, "missense_variant"),
    (761139, "C", "T", 0.62, "missense_variant"),
    (761155, "C", "T", 0.97, "missense_variant"),
    (761161, "T", "C", 0.12, "missense_variant"),
    (762368, "G", "A", 1.00, "synonymous_variant"),
    (762917, "C", "G", 0.24, "missense_variant"),
]


def path(name):
    return os.path.join(HERE, name)


def reference(rng):
    # A GC content near the 65% this genome sits at.
    bases = []
    for _ in range(START, END + 1):
        bases.append(rng.choice("GGCC" * 13 + "AATT" * 7))
    # The reference base at each call is the one the call names.
    for position, ref, _, _, _ in CALLS:
        bases[position - START] = ref
    return "".join(bases)


def write_reference(seq):
    with open(path("ref.fa"), "w") as out:
        # A slice written the way samtools faidx writes one, so a reader knows
        # where it starts.
        out.write(f">{SEQ}:{START}-{END}\n")
        for at in range(0, len(seq), 60):
            out.write(seq[at:at + 60] + "\n")


def write_genes():
    with open(path("genes.gff3"), "w") as out:
        out.write("##gff-version 3\n")
        for name, start, end, strand in GENES:
            out.write(f"{SEQ}\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={name};Name={name}\n")


def write_calls():
    plain = path("calls.vcf")
    with open(plain, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={SEQ},length={LENGTH}>\n")
        out.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        out.write('##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotation">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for position, ref, alt, af, kind in CALLS:
            out.write(
                f"{SEQ}\t{position}\t.\t{ref}\t{alt}\t60\tPASS\t"
                f"AF={af:.2f};ANN={alt}|{kind}|MODERATE|rpoB\n"
            )
    subprocess.run(["bgzip", "-f", plain], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", plain + ".gz"], check=True)


def write_reads(rng, seq):
    sam = path("reads.sam")
    with open(sam, "w") as out:
        out.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{SEQ}\tLN:{LENGTH}\n")
        count = 0
        for start in range(START, END - 100 + 2, 4):
            # Fewer reads over a stretch the sample has lost, so the depth
            # drops there as a deletion leaves it.
            if 761880 <= start + 50 <= 762040 and rng.random() < 0.9:
                continue
            if rng.random() < 0.2:
                continue
            read = list(seq[start - START:start - START + 100])
            for position, _, alt, af, _ in CALLS:
                offset = position - start
                if 0 <= offset < 100 and rng.random() < af:
                    read[offset] = alt
            count += 1
            flag = 16 if rng.random() < 0.5 else 0
            out.write(
                f"read{count:05d}\t{flag}\t{SEQ}\t{start}\t60\t100M\t*\t0\t0\t"
                f"{''.join(read)}\t{'I' * 100}\n"
            )
    subprocess.run(["samtools", "sort", "-o", path("reads.bam"), sam], check=True)
    subprocess.run(["samtools", "index", path("reads.bam")], check=True)
    os.remove(sam)


def write_scan(rng):
    # One chromosome as PLINK writes it, called 1, with a peak over rpoB.
    with open(path("gwas.assoc"), "w") as out:
        out.write(" CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P        OR\n")
        for number, bp in enumerate(range(2000, LENGTH, 2200)):
            distance = abs(bp - 761155)
            # The signal fades with distance from the causal site, and each
            # marker keeps only part of it, as markers in partial linkage do.
            strength = max(0.0, 11.0 - distance / 9000.0) * (0.25 + 0.75 * rng.random())
            p = min(1.0, rng.random() * 10 ** (-strength))
            p = max(p, 1e-300)
            out.write(
                f"   1   snp{number:05d} {bp:10d}    T   0.3000   0.2500    C        5.000 {p:12.4g}     1.285\n"
            )


def write_samples():
    # Forty samples in four lineages, a tree of them with branch lengths and
    # support, a second tree that disagrees in two places, a sample sheet in
    # an order that is not the tree's, and an alignment evolved down the tree.
    # Its own seed, so these files do not move when the ones above change.
    rng = random.Random(20260925)
    lineages = {"L1": 6, "L2": 10, "L3": 8, "L4": 16}
    countries = {"L1": ["Vietnam", "India"], "L2": ["Vietnam", "China", "Peru"],
                 "L3": ["India", "Kenya"], "L4": ["Spain", "Peru", "Kenya", "Portugal"]}
    samples = []
    number = 0
    for lineage, size in lineages.items():
        for _ in range(size):
            number += 1
            samples.append({"sample": f"S{number:02d}", "lineage": lineage,
                            "country": rng.choice(countries[lineage]),
                            "year": rng.randint(2008, 2023),
                            "host": rng.choice(["human"] * 5 + ["cattle"])})

    def build(tips):
        nodes = list(tips)
        rng.shuffle(nodes)
        while len(nodes) > 1:
            a = nodes.pop(rng.randrange(len(nodes)))
            b = nodes.pop(rng.randrange(len(nodes)))
            nodes.append((a, b))
        return nodes[0]

    clades = {lineage: build([s["sample"] for s in samples if s["lineage"] == lineage])
              for lineage in lineages}
    tree = (((clades["L1"], clades["L2"]), clades["L3"]), clades["L4"])

    def newick(node, depth=0):
        if isinstance(node, str):
            return f"{node}:{round(rng.uniform(0.002, 0.02), 4)}"
        left, right = node
        inner = f"({newick(left, depth + 1)},{newick(right, depth + 1)})"
        length = round(rng.uniform(0.001, 0.015), 4)
        label = str(rng.randint(55, 100))
        return f"{inner}{label}:{length}" if depth else inner

    with open(path("tree.nwk"), "w") as out:
        out.write(newick(tree) + ";\n")
    moved = [s["sample"] for s in samples if s["lineage"] == "L2"][:2]
    rest = [s["sample"] for s in samples if s["lineage"] == "L2"][2:]
    first = build([s["sample"] for s in samples if s["lineage"] == "L1"] + moved)
    second = build(rest)
    with open(path("tree2.nwk"), "w") as out:
        out.write(newick(((first, second), (clades["L3"], clades["L4"]))) + ";\n")
    sheet = samples[:]
    rng.shuffle(sheet)
    with open(path("samples.tsv"), "w") as out:
        out.write("sample\tlineage\tcountry\tyear\thost\n")
        for s in sheet:
            out.write(f"{s['sample']}\t{s['lineage']}\t{s['country']}\t{s['year']}\t{s['host']}\n")
    root = [rng.choice("ACGT") for _ in range(300)]
    rows = {}

    def evolve(node, seq):
        seq = seq[:]
        for _ in range(rng.randint(1, 6)):
            seq[rng.randrange(len(seq))] = rng.choice("ACGT")
        if isinstance(node, str):
            rows[node] = "".join(seq)
        else:
            evolve(node[0], seq)
            evolve(node[1], seq)

    evolve(tree, root)
    with open(path("aln.fasta"), "w") as out:
        out.write("".join(f">{name}\n{rows[name]}\n" for name in sorted(rows)))


def write_assemblies(rng):
    # Two assemblies of one chromosome, the second with a stretch turned round
    # and a stretch moved, aligned as minimap2 writes it.
    qlen, tlen = 600_000, 620_000
    pieces = [  # query start, query end, strand, target start
        (0, 150_000, "+", 0),
        (150_000, 260_000, "-", 150_000),
        (260_000, 380_000, "+", 400_000),
        (380_000, 470_000, "+", 280_000),
        (470_000, 600_000, "+", 520_000),
    ]
    with open(path("assemblies.paf"), "w") as out:
        for qstart, qend, strand, tstart in pieces:
            # Each piece in stretches, as an aligner breaks a long one at its
            # gaps, a little less than all of it matching.
            at = qstart
            while at < qend:
                step = min(qend - at, rng.randint(15_000, 40_000))
                span = step - rng.randint(0, 400)
                offset = at - qstart
                if strand == "+":
                    t0 = tstart + offset
                else:
                    t0 = tstart + (qend - qstart) - offset - span
                matches = int(span * rng.uniform(0.96, 0.995))
                out.write(
                    f"asm1_chr1\t{qlen}\t{at}\t{at + span}\t{strand}\tasm2_chr1\t{tlen}\t"
                    f"{t0}\t{t0 + span}\t{matches}\t{span}\t60\ttp:A:P\n"
                )
                at += step


def write_windows(rng):
    # The depth of two samples along a whole chromosome in windows of 10 kb:
    # the first has lost a stretch and carries another twice, the second is
    # whole and was sequenced deeper.
    samples = {
        "sampleA.bedgraph": (48, [(1_450_000, 1_510_000, 0.02), (3_120_000, 3_180_000, 2.0)]),
        "sampleB.bedgraph": (70, []),
    }
    for name, (mean, changes) in samples.items():
        with open(path(name), "w") as out:
            level = 0.0
            for start in range(0, LENGTH, 10_000):
                end = min(start + 10_000, LENGTH)
                level = 0.8 * level + rng.gauss(0, 1.5)
                value = mean + level
                for low, high, factor in changes:
                    if low <= start < high:
                        value *= factor
                out.write(f"{SEQ}\t{start}\t{end}\t{max(0.0, value):.1f}\n")


def write_zip():
    # Every file a page draws, in one download. Dated the same every time, so
    # the archive is the same file when its contents are.
    import zipfile
    names = ["reads.bam", "reads.bam.bai", "genes.gff3", "calls.vcf.gz",
             "calls.vcf.gz.tbi", "ref.fa", "gwas.assoc", "tree.nwk", "tree2.nwk",
             "samples.tsv", "aln.fasta", "assemblies.paf", "sampleA.bedgraph",
             "sampleB.bedgraph"]
    with zipfile.ZipFile(path("examples.zip"), "w", zipfile.ZIP_DEFLATED) as out:
        for name in names:
            info = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            with open(path(name), "rb") as held:
                out.writestr(info, held.read())


def main():
    rng = random.Random(20260926)
    seq = reference(rng)
    write_reference(seq)
    write_genes()
    write_calls()
    write_reads(rng, seq)
    write_scan(rng)
    write_samples()
    write_assemblies(rng)
    write_windows(rng)
    write_zip()


if __name__ == "__main__":
    main()
