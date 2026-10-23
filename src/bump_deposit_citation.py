#!/usr/bin/env python3
"""Rewrite every citation of the archival deposit, from one command.

The deposit this version of the work is archived as is one fact written out by
hand in seven places across six files. Bumping it by hand after each release is
how a document comes to point at an older deposit than its neighbours. This
rewrites all seven, and src/check_manuscript.py asserts afterwards that they
agree with docs/archive_snapshot.tsv.

    python3 src/bump_deposit_citation.py v7.0 10.5281/zenodo.NNNNNNNN
    make archive-snapshot TAG=v7.0 DOI=10.5281/zenodo.NNNNNNNN

Run in that order: this rewrites the citations, the snapshot records what the
deposit contains, and the gate then checks the two against each other.
"""

import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DOCS = ("manuscript/manuscript_v6.0.md",
        "manuscript/supplementary_v6.0.md",
        "manuscript/title_page.md",
        "manuscript/cover_letter.md",
        "README.md")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    tag, doi = sys.argv[1], sys.argv[2]
    if not re.fullmatch(r"v\d+\.\d+", tag):
        sys.exit("not a version tag: %r" % tag)
    m = re.fullmatch(r"10\.5281/zenodo\.(\d+)", doi)
    if not m:
        sys.exit("not a Zenodo version DOI: %r" % doi)
    num = m.group(1)

    changed = 0
    for rel in DOCS:
        p = os.path.join(ROOT, rel)
        if not os.path.isfile(p):
            continue
        body = io.open(p, encoding="utf-8").read()
        out = []
        for line in body.split("\n"):
            # The concept DOI names the series, not a version, and is left alone
            # wherever the line says so.
            if "concept" not in line.lower():
                line = re.sub(r"zenodo\.\d+", "zenodo." + num, line)
            line = re.sub(r"(at tag )v\d+\.\d+", r"\g<1>" + tag, line)
            line = re.sub(r"(\(version )v\d+\.\d+", r"\g<1>" + tag, line)
            out.append(line)
        new = "\n".join(out)
        if new != body:
            io.open(p, "w", encoding="utf-8").write(new)
            changed += 1
            print("rewrote " + rel)

    # CITATION.cff keeps superseded identifiers by design. The current deposit is
    # inserted immediately after the concept entry rather than replacing the list.
    cff = os.path.join(ROOT, "CITATION.cff")
    if os.path.isfile(cff):
        c = io.open(cff, encoding="utf-8").read()
        want = ("  - type: doi\n    value: %s\n"
                "    description: Archived snapshot of version %s\n" % (doi, tag))
        if want in c:
            print("CITATION.cff already lists " + tag)
        else:
            anchor = re.search(
                r"  - type: doi\n    value: [^\n]+\n"
                r"    description: Concept DOI[^\n]*\n", c)
            if not anchor:
                sys.exit("CITATION.cff: no concept DOI entry to insert after")
            c = c[:anchor.end()] + want + c[anchor.end():]
            io.open(cff, "w", encoding="utf-8").write(c)
            changed += 1
            print("rewrote CITATION.cff")

    print("%d files rewritten; run make archive-snapshot TAG=%s DOI=%s next"
          % (changed, tag, doi))


if __name__ == "__main__":
    main()
