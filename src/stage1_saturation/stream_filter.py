#!/usr/bin/env python3
"""Stream a gzipped GWAS summary statistics file from stdin and keep only rows
matching the committed union instrument set. Reproduces src/stage1_saturation/extract.py
stream()/byrs()/bykey()/fg() exactly, but reads the union set from the committed
union_snps.tsv instead of rebuilding it from EUR.bim, and reads the compressed
stream from stdin so nothing large is written to disk.

usage: zcat file.gz | stream_filter.py <union_snps.tsv> <mode> <out.tsv> [args]
  mode byrs  <colname>              keep rows whose <colname> is a union rsid
  mode bykey <chrcol> <poscol>      keep rows whose chr:pos maps to a union rsid
  mode fg                           keep rows whose comma separated rsids hit the union
"""
import sys

union, mode, out = sys.argv[1], sys.argv[2], sys.argv[3]

rs = set()
keys = {}
with open(union) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    i_chr, i_rsid, i_pos = h.index("chr"), h.index("rsid"), h.index("pos")
    for line in fh:
        x = line.rstrip("\n").split("\t")
        rs.add(x[i_rsid])
        keys[x[i_chr] + ":" + x[i_pos]] = x[i_rsid]

fh = sys.stdin
o = open(out, "w")
head = fh.readline().rstrip("\r\n")
H = head.split("\t")
o.write(head + "\tSNP\n")
kept = 0
seen = 0

if mode == "byrs":
    j = H.index(sys.argv[4])
    for line in fh:
        seen += 1
        x = line.rstrip("\r\n").split("\t")
        if x[j] in rs:
            o.write(line.rstrip("\r\n") + "\t" + x[j] + "\n")
            kept += 1
elif mode == "bykey":
    jc, jp = H.index(sys.argv[4]), H.index(sys.argv[5])
    for line in fh:
        seen += 1
        x = line.rstrip("\r\n").split("\t")
        s = keys.get(x[jc] + ":" + x[jp])
        if s:
            o.write(line.rstrip("\r\n") + "\t" + s + "\n")
            kept += 1
elif mode == "fg":
    j = H.index("rsids")
    for line in fh:
        seen += 1
        x = line.rstrip("\r\n").split("\t")
        for v in x[j].split(","):
            if v in rs:
                o.write(line.rstrip("\r\n") + "\t" + v + "\n")
                kept += 1
                break
else:
    sys.exit("unknown mode " + mode)

o.close()
sys.stderr.write("%s: scanned %d rows, kept %d\n" % (out, seen, kept))
