#!/bin/sh
# Writes the .hic files the reader beside this folder is tested on, with
# hictk, and then reads windows of them back with hictk dump. The reader's
# tests compare what it reads with those dumps, so it is held to what hictk
# says a file holds, not to what this folder's author thought it should hold.
#
# The tool is hictk 2.2.0 from bioconda, which writes version 9, the version
# the reader reads. Run from anywhere:
#
#     sh src/read/fixtures/hic/make.sh
#
# The contacts are written by a rule, not drawn at random, so the file is the
# same wherever it is made. They are sparse on purpose: a few hundred cells
# of a chromosome of 1,235 bins, which is past the thousand bins at which
# hictk cuts a map into more than one column of blocks, and contacts more
# than 874 bins apart, which hictk files in the second band of blocks from
# the diagonal. Every block a large map has to be found by is there, in a
# file of a few kilobytes.
set -eu
cd "$(dirname "$0")"

printf 'chr1\t1234567\nchr2\t345678\n' > contacts.sizes

awk -v OFS='\t' '
function cell(c, size, i, j, value) {
    print c, i * 1000, ends(i, size), c, j * 1000, ends(j, size), value
}
function ends(bin, size) {
    return (bin + 1) * 1000 > size ? size : (bin + 1) * 1000
}
BEGIN {
    for (i = 0; i < 1235; i++) {
        if (i % 4 == 0) cell("chr1", 1234567, i, i, 1 + (i * 37) % 50)
        if (i % 6 == 1 && i + 1 < 1235) cell("chr1", 1234567, i, i + 1, (1 + (i * 11) % 20) + 0.5)
        if (i % 97 == 5) for (j = i + 880; j < 1235; j += 89) cell("chr1", 1234567, i, j, 1 + (i + j) % 7)
    }
    for (i = 0; i < 346; i++) {
        if (i % 3 == 0) cell("chr2", 345678, i, i, 2 + i % 9)
        if (i % 5 == 2 && i + 2 < 346) cell("chr2", 345678, i, i + 2, 1)
    }
    for (i = 0; i < 20; i++) print "chr1", i * 1000, (i + 1) * 1000, "chr2", i * 1000, (i + 1) * 1000, 3
}' > contacts.bg2

hictk load -f bg2 -c contacts.sizes -b 1000 --count-as-float --force -v 1 \
    contacts.bg2 base.hic
hictk zoomify --resolutions 1000 2000 5000 10000 50000 250000 --force -v 1 \
    base.hic contacts.hic
rm base.hic contacts.bg2 contacts.sizes

# Each window as `# RESOLUTION SEQUENCE START END`, 0-based and half-open as
# hictk takes a range, then the cells hictk dump prints for it.
for window in '1000 chr1 300000 420000' '1000 chr1 0 1234567' \
    '2000 chr1 1100000 1234567' '5000 chr1 612345 987654' \
    '50000 chr1 0 1234567' '250000 chr1 0 1234567' '10000 chr2 0 345678' \
    '250000 chr2 0 345678' '1000 chr2 100000 100001'; do
    set -- $window
    echo "# $window"
    hictk dump --join --resolution "$1" -r "$2:$3-$4" contacts.hic
done > contacts.dump

# A second file, of one square of contacts far from the diagonal: the first
# ten bins of a sequence against bins 40 to 49, every cell a different count.
# hictk writes it as a dense block whose columns start at bin 0 and whose
# rows start at bin 40, which no block of contacts.hic is: a dense block on
# the diagonal is its own mirror, and a cell read with its column and its row
# swapped lands where its mirror was. Here it lands on another cell's count.
printf 'chr1\t56789\n' > square.sizes
awk -v OFS='\t' 'BEGIN {
    for (i = 0; i < 10; i++)
        for (j = 40; j < 50; j++)
            print "chr1", i * 1000, (i + 1) * 1000, "chr1", j * 1000, \
                (j + 1) * 1000, 1 + 10 * i + j - 40
}' > square.bg2
hictk load -f bg2 -c square.sizes -b 1000 --count-as-float --force -v 1 \
    square.bg2 square.hic
rm square.bg2 square.sizes

# The whole sequence, and a window whose edges cut through the square.
for window in '1000 chr1 0 56789' '1000 chr1 5500 44500'; do
    set -- $window
    echo "# $window"
    hictk dump --join --resolution "$1" -r "$2:$3-$4" square.hic
done > square.dump
