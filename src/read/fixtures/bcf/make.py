#!/usr/bin/env python3
"""Writes the BCF files the BCF reader and the command line are tested on, and
what bcftools prints for each, from the VCF text beside them.

A reader's test compares what it writes for a BCF with what bcftools prints
for the same file, so it is held to what htslib says the file holds, not to
what this script thought it should hold. The tools are htslib 1.24's
(bcftools and bgzip). Run from anywhere:

    python3 src/read/fixtures/bcf/make.py

- tiny.vcf holds a value of every type a BCF stores and every way a value is
  missing; tiny.bcf is it as `bcftools view -Ob` writes it, with the CSI
  `bcftools index` writes, tiny.u.bcf as `-Ou` writes it, BGZF of blocks
  stored rather than compressed, and tiny.raw.bcf the bare stream inside.
- floats.bcf holds floats of every size, a few chosen for where htslib's
  rounding and printf's part ways and the rest drawn from a seed, less the
  ones printf prints differently on macOS and on Linux.
- phased.vcf is VCF 4.4, whose first allele may say how it is phased.
- cohort.bcf holds the calls of ../indexed/cohort.vcf.gz, and sv.bcf the
  structural calls of sv.vcf, breakends whose mates lie far off among them,
  each in blocks of a few kilobytes cut anywhere, a record across two blocks included, as make.py there cuts
  its files, so a window is read out of the middle of a file of many blocks
  as it is out of a large one.
"""

import decimal
import gzip
import os
import random
import struct
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
INDEXED = os.path.join(HERE, "..", "indexed")
# The empty block bgzip ends a file with, left off every piece but the last.
EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def run(*command, out=None):
    result = subprocess.run(command, cwd=HERE, check=True, capture_output=True)
    if out is not None:
        with open(os.path.join(HERE, out), "wb") as file:
            file.write(result.stdout)


def printed(name):
    """What bcftools prints for `name`, every field of every record."""
    run("bcftools", "view", "--no-version", name, out=name + ".vcf")


def blocks(raw, name, piece):
    """The bare BCF stream `raw` in BGZF blocks of about `piece` bytes, cut
    anywhere."""
    out = bytearray()
    for at in range(0, len(raw), piece):
        block = subprocess.run(
            ["bgzip", "-c"], input=raw[at : at + piece], capture_output=True, check=True
        ).stdout
        assert block.endswith(EOF)
        out += block[: -len(EOF)]
    out += EOF
    with open(os.path.join(HERE, name), "wb") as file:
        file.write(out)


def raw_of(vcf):
    """The bare BCF stream bcftools writes for `vcf`."""
    stored = subprocess.run(
        ["bcftools", "view", "--no-version", "-Ou", vcf],
        cwd=HERE,
        check=True,
        capture_output=True,
    ).stdout
    return gzip.decompress(stored)


def tiny():
    run("bcftools", "view", "--no-version", "-Ob", "-o", "tiny.bcf", "tiny.vcf")
    run("bcftools", "index", "-f", "tiny.bcf")
    run("bcftools", "view", "--no-version", "-Ou", "-o", "tiny.u.bcf", "tiny.vcf")
    with open(os.path.join(HERE, "tiny.raw.bcf"), "wb") as file:
        file.write(raw_of("tiny.vcf"))
    printed("tiny.bcf")
    run("bcftools", "view", "--no-version", "-G", "tiny.bcf", out="tiny.bcf.sites.vcf")
    run(
        "bcftools",
        "annotate",
        "--no-version",
        "-x",
        "^FORMAT/GT",
        "tiny.bcf",
        out="tiny.bcf.gt.vcf",
    )


def tie(text):
    """Whether `text`, stored as a 32-bit float, is printed by printf, which
    htslib hands values outside 0.0001 to 999,999, as a half to round: macOS's
    libc keeps a trailing nought after such a rounding where glibc drops it,
    so the text bcftools prints for one depends on the system it ran on."""
    stored = struct.unpack("<f", struct.pack("<f", float(text)))[0]
    if 0.0001 <= abs(stored) <= 999999 or stored == 0:
        return False
    digits = decimal.Decimal(stored).normalize().as_tuple().digits
    return len(digits) == 7 and digits[-1] == 5


def floats():
    chosen = [
        "0.0001", "0.00012345", "0.000099999", "0.00009999995", "1e-05", "1.5e-07",
        "0.333333", "0.1", "0.25", "0.5", "2.5", "0.7", "0.07", "0.007", "0.0007",
        "999999", "999999.4", "999999.6", "1000000", "1234565", "1234575", "123456.5",
        "99999.95", "12345.65", "1.0000005", "7", "-0.5", "-0", "0", "3.4e38",
        "1.17549e-38", "16777216", "0.1234565", "65.4321", "8.88888",
    ]
    rng = random.Random(20261002)
    drawn = []
    while len(drawn) < 600:
        exponent = rng.uniform(-9, 9)
        value = f"{rng.choice([1, -1]) * 10**exponent:.9g}"
        if not tie(value):
            drawn.append(value)
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=100000>",
        '##INFO=<ID=F,Number=.,Type=Float,Description="Floats">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    values = chosen + drawn
    for at in range(0, len(values), 20):
        group = values[at : at + 20]
        quality = group[0].lstrip("-")
        lines.append(f"chr1\t{at + 1}\t.\tA\tG\t{quality}\t.\tF={','.join(group)}")
    with open(os.path.join(HERE, "floats.vcf"), "w") as file:
        file.write("\n".join(lines) + "\n")
    run("bcftools", "view", "--no-version", "-Ob", "-o", "floats.bcf", "floats.vcf")
    printed("floats.bcf")


def phased():
    run("bcftools", "view", "--no-version", "-Ob", "-o", "phased.bcf", "phased.vcf")
    printed("phased.bcf")


def cut(vcf, name, piece):
    blocks(raw_of(vcf), name, piece)
    run("bcftools", "index", "-f", name)
    printed(name)


def main():
    for tool in ("bcftools", "bgzip"):
        version = subprocess.run([tool, "--version"], capture_output=True, text=True).stdout
        print(version.splitlines()[0])
    tiny()
    floats()
    phased()
    cut(os.path.join(INDEXED, "cohort.vcf.gz"), "cohort.bcf", 3000)
    cut("sv.vcf", "sv.bcf", 200)


if __name__ == "__main__":
    main()
