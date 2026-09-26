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


def write_trait_scan():
    # A scan across a whole genome of twelve chromosomes, as PLINK writes one,
    # for the figure of every chromosome at once. It keeps a seed of its own,
    # so the files above stay as they are when this changes. Two peaks, on 3
    # and on 9, cross the genome-wide line.
    rng = random.Random(20261001)
    lengths = [43, 36, 36, 35, 30, 31, 30, 28, 23, 23, 29, 27]  # megabases
    peaks = {3: 18_500_000, 9: 11_200_000}
    with open(path("trait.assoc"), "w") as out:
        out.write(" CHR          SNP         BP   A1      F_A      F_U   A2        CHISQ            P        OR\n")
        number = 0
        for chrom, megabases in enumerate(lengths, start=1):
            for bp in range(40_000, megabases * 1_000_000, 150_000):
                strength = 0.0
                if chrom in peaks:
                    strength = max(0.0, 10.0 - abs(bp - peaks[chrom]) / 120_000.0)
                strength *= 0.25 + 0.75 * rng.random()
                p = max(min(1.0, rng.random() * 10 ** (-strength)), 1e-300)
                out.write(
                    f"{chrom:4d}   snp{number:05d} {bp:10d}    T   0.3000   0.2500    C        5.000 {p:12.4g}     1.285\n"
                )
                number += 1


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


def write_lineages(rng):
    # Weekly counts of four lineages out of every genome sequenced that week:
    # A fades, B.1 rises and is overtaken by B.2, and C stays rare. Each
    # lineage's share follows its own growth rate, and the counts are drawn
    # from it, so the frequencies wander the way sampled ones do.
    import math
    growth = {"A": (3.0, -0.18), "B.1": (0.4, 0.13), "B.2": (-4.5, 0.36), "C": (0.3, 0.0)}
    with open(path("lineages.tsv"), "w") as out:
        out.write("week\tlineage\tcount\ttotal\n")
        for week in range(1, 31):
            total = 60 + int(40 * (1 + math.sin(week / 4))) + rng.randrange(0, 20)
            scores = {name: a + b * week for name, (a, b) in growth.items()}
            top = max(scores.values())
            weights = {name: math.exp(score - top) for name, score in scores.items()}
            whole = sum(weights.values())
            left = total
            names = list(growth)
            for index, name in enumerate(names):
                if index == len(names) - 1:
                    count = left
                else:
                    share = weights[name] / whole
                    count = min(left, sum(1 for _ in range(total) if rng.random() < share))
                left -= count
                out.write(f"{week}\t{name}\t{count}\t{total}\n")


def write_reproduction(rng):
    # The reproductive number over the same weeks, with its 95% interval,
    # as EpiEstim or a phylodynamic model writes one: above one while B.2
    # spreads, below it between the waves. The interval narrows where more
    # cases were seen.
    import math
    with open(path("reproduction.tsv"), "w") as out:
        out.write("week\tmean\tlower\tupper\n")
        for week in range(3, 31):
            mean = 1.0 + 0.32 * math.sin((week - 2) / 4.2) + rng.gauss(0, 0.03)
            width = 0.12 + 0.25 / math.sqrt(week)
            out.write(f"{week}\t{mean:.3f}\t{mean - width:.3f}\t{mean + width * 1.2:.3f}\n")


def write_selection(rng):
    # A test of selection at each of 300 codons, as HyPhy's FEL writes it:
    # most sites under purifying selection, beta below alpha, and a handful
    # in two stretches where beta is well above it and the test says so.
    with open(path("fel.csv"), "w") as out:
        out.write("site,alpha,beta,p-value\n")
        hot = set(range(58, 66)) | set(range(181, 187))
        for site in range(1, 301):
            alpha = rng.uniform(0.3, 2.0)
            if site in hot and rng.random() < 0.7:
                beta = alpha * rng.uniform(3.0, 9.0)
                p = rng.uniform(0.0005, 0.04)
            elif rng.random() < 0.08:
                beta = alpha * rng.uniform(0.9, 2.0)
                p = rng.uniform(0.1, 0.9)
            else:
                beta = alpha * rng.uniform(0.02, 0.6)
                p = rng.uniform(0.2, 1.0)
            out.write(f"{site},{alpha:.3f},{beta:.3f},{p:.4f}\n")


def write_signal(rng):
    # Two nanopore reads as slow5tools view writes them: the current steps
    # from level to level as the strand ratchets through the pore, each step
    # a few samples to a few dozen long, with noise on top. Stored raw, and
    # put into picoamperes with each read's digitisation, offset and range.
    #
    # And the basecaller's record of each, as Dorado writes it with
    # --emit-moves: a base for each step, its move table in strides of five
    # samples with a 1 where a step starts. The bases come from a generator
    # of their own, so the signals are what they were before there were any.
    digitisation, offset, span = 8192.0, 6.0, 1467.61
    stride = 5
    called = random.Random(5)
    records = []
    with open(path("reads.slow5"), "w") as out:
        out.write("#slow5_version\t0.2.0\n#num_read_groups\t1\n")
        out.write("@run_id\texample\n")
        out.write("#char*\tuint32_t\tdouble\tdouble\tdouble\tdouble\tuint64_t\tint16_t*\n")
        out.write("#read_id\tread_group\tdigitisation\toffset\trange\tsampling_rate"
                  "\tlen_raw_signal\traw_signal\n")
        for read, samples in (("read_1", 2400), ("read_2", 1800)):
            raw = []
            starts = []
            level = rng.uniform(80, 100)
            while len(raw) < samples:
                level = min(125.0, max(65.0, level + rng.gauss(0, 12)))
                starts.append(len(raw))
                for _ in range(rng.randint(4, 40)):
                    current = level + rng.gauss(0, 1.8)
                    raw.append(round(current * digitisation / span - offset))
            raw = raw[:samples]
            out.write(f"{read}\t0\t{digitisation:.0f}\t{offset:.0f}\t{span}\t4000"
                      f"\t{samples}\t{','.join(str(value) for value in raw)}\n")
            begun = {start // stride for start in starts if start < samples}
            moves = [1 if block in begun else 0 for block in range(samples // stride)]
            moves[0] = 1
            seq = "".join(called.choice("ACGT") for _ in range(sum(moves)))
            records.append(f"{read}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t*\t"
                           f"mv:B:c,{stride},{','.join(map(str, moves))}\tts:i:0\tns:i:{samples}\n")
    with open(path("moves.sam"), "w") as out:
        out.write("@HD\tVN:1.6\tSO:unknown\n@PG\tID:basecaller\tPN:dorado\n")
        out.writelines(records)


def write_depths(rng):
    # The depth of all 40 samples of the tree in windows of 100 kb along the
    # whole chromosome, as bedtools unionbedg writes it. Each sample was
    # sequenced to its own depth; one clade has lost a 100 kb stretch and
    # another carries a stretch twice.
    import re
    tips = re.findall(r"(S\d\d):", open(path("tree.nwk")).read())
    lost = {"S08", "S10", "S11", "S12", "S13", "S14", "S15", "S16"}
    doubled = {"S18", "S23", "S24"}
    mean = {tip: rng.uniform(45, 110) for tip in tips}
    with open(path("depths.tsv"), "w") as out:
        out.write("chrom\tstart\tend\t" + "\t".join(tips) + "\n")
        for start in range(0, LENGTH, 100_000):
            end = min(start + 100_000, LENGTH)
            values = []
            for tip in tips:
                depth = mean[tip] * rng.uniform(0.9, 1.1)
                if tip in lost and 1_400_000 <= start < 1_500_000:
                    depth = 0.0
                if tip in doubled and 3_100_000 <= start < 3_200_000:
                    depth *= 2
                values.append(f"{depth:.1f}")
            out.write(f"{SEQ}\t{start}\t{end}\t" + "\t".join(values) + "\n")


def write_linkage(rng):
    # Linkage between 36 variants across rpoB, as PLINK's --r2 writes it with
    # --ld-window-r2 0: three blocks inherited together, strong inside each
    # and weak between them, fading with distance.
    import math
    sites = sorted(rng.sample(range(759_900, 763_300), 36))
    blocks = [760_900, 762_100]  # where one block ends and the next begins
    block = [sum(site >= edge for edge in blocks) for site in sites]
    with open(path("linkage.ld"), "w") as out:
        out.write(" CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2 \n")
        for i, a in enumerate(sites):
            for j in range(i + 1, len(sites)):
                b = sites[j]
                if block[i] == block[j]:
                    r2 = 0.95 * math.exp(-(b - a) / 4000) * rng.uniform(0.75, 1.0)
                else:
                    r2 = 0.25 * math.exp(-(b - a) / 1500) * rng.uniform(0.0, 1.0)
                out.write(f" {SEQ} {a:>12} {f'v{i + 1}':>12} {SEQ} {b:>12} {f'v{j + 1}':>12} "
                          f"{r2:>12.4f} \n")


def write_epistasis():
    # A few pairs of the calls that change together, with how strongly: a
    # table of your own, headed by what its columns are.
    pairs = [
        (761110, 761155, 0.82), (761139, 761161, 0.35), (760314, 762368, 0.64),
        (761155, 762917, 0.91), (761110, 762917, 0.28),
    ]
    with open(path("epistasis.tsv"), "w") as out:
        out.write("pos1\tpos2\tscore\n")
        for a, b, score in pairs:
            out.write(f"{a}\t{b}\t{score}\n")


def write_lead_linkage(rng):
    # The linkage of every marker within 100 kb of the scan's strongest with
    # it, as PLINK's --r2 --ld-snp writes it, and a genetic map of the same
    # stretch every 5 kb, as HapMap writes one. A marker is linked in
    # proportion to how much of the signal it carries, as markers in partial
    # linkage are.
    import math
    rows = []
    with open(path("gwas.assoc")) as held:
        next(held)
        for line in held:
            fields = line.split()
            rows.append((fields[1], int(fields[2]), float(fields[8])))
    lead_name, lead, lead_p = min(rows, key=lambda row: row[2])
    top = -math.log10(lead_p)
    with open(path("lead.ld"), "w") as out:
        out.write(" CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2 \n")
        for name, bp, p in rows:
            if name == lead_name or abs(bp - lead) > 100_000:
                continue
            carried = min(1.0, max(0.0, -math.log10(p) / top))
            r2 = carried ** 0.7 * rng.uniform(0.85, 1.0)
            out.write(f"     1 {lead:>12} {lead_name:>12}      1 {bp:>12} {name:>12} {r2:>12.4f} \n")
    with open(path("genetic_map.txt"), "w") as out:
        out.write("Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n")
        distance = 0.0
        for start in range(lead - 150_000, lead + 150_000, 5_000):
            rate = 0.4 + rng.uniform(0, 0.6)
            for hotspot, height in ((lead - 62_000, 38.0), (lead + 47_000, 55.0)):
                rate += height * math.exp(-((start - hotspot) / 6_000) ** 2)
            rate = round(rate, 2)
            out.write(f"1\t{start + 1}\t{rate:.2f}\t{distance:.6f}\n")
            distance += rate * 5_000 / 1e6


def write_zip():
    # Every file a page draws, in one download. Dated the same every time, so
    # the archive is the same file when its contents are.
    import zipfile
    names = ["reads.bam", "reads.bam.bai", "genes.gff3", "calls.vcf.gz",
             "calls.vcf.gz.tbi", "ref.fa", "gwas.assoc", "tree.nwk", "tree2.nwk",
             "samples.tsv", "aln.fasta", "assemblies.paf", "sampleA.bedgraph",
             "sampleB.bedgraph", "lineages.tsv", "reproduction.tsv", "fel.csv",
             "reads.slow5", "moves.sam", "depths.tsv", "linkage.ld", "epistasis.tsv",
             "lead.ld", "genetic_map.txt", "trait.assoc"]
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
    write_lineages(rng)
    write_reproduction(rng)
    write_selection(rng)
    write_signal(rng)
    write_depths(rng)
    write_linkage(rng)
    write_epistasis()
    write_lead_linkage(rng)
    write_trait_scan()
    write_zip()


if __name__ == "__main__":
    main()
