"""Write docs/archive_snapshot.tsv: what the cited archival deposit contains.

A deposit is cut from a tag, and the commit that updates the manuscript to cite
that tag and its DOI necessarily comes after it. That one commit is the only
permitted difference between the archive and the submitted manuscript, so the
hashes recorded here are computed after masking the archive citation itself:
every DOI and every version tag is replaced by a token before hashing. A
citation update therefore leaves the hash unchanged, while any other edit to a
document the archive is supposed to contain changes it.

usage: python3 src/write_archive_snapshot.py v6.3 10.5281/zenodo.23326289
"""

import hashlib
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COVERED = ["manuscript/manuscript_v6.0.md",
           "manuscript/supplementary_v6.0.md",
           "manuscript/STROBE-MR_checklist.md"]


def masked(path):
    body = open(path, encoding="utf-8").read()
    body = re.sub(r"10\.5281/zenodo\.\d+", "<DOI>", body)
    body = re.sub(r"\bv\d+\.\d+(\.\d+)?\b", "<TAG>", body)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    tag, doi = sys.argv[1], sys.argv[2]
    if not re.fullmatch(r"v\d+\.\d+(\.\d+)?", tag):
        sys.exit("the tag must look like a version tag, got %r" % tag)
    if not re.fullmatch(r"10\.5281/zenodo\.\d+", doi):
        sys.exit("the DOI must be a Zenodo version DOI, got %r" % doi)
    out = [
        "# What the archival deposit the manuscript cites contains.",
        "# Regenerate with: make archive-snapshot TAG=vX.Y DOI=10.5281/zenodo.N",
        "# Hashes are of the file with every DOI and version tag masked, so the",
        "# commit that updates the citation after the tag is cut cannot register",
        "# as drift while any other change does.",
        "tag\t%s" % tag,
        "version_doi\t%s" % doi,
    ]
    for rel in COVERED:
        out.append("file\t%s\t%s" % (rel, masked(os.path.join(ROOT, rel))))
    path = os.path.join(ROOT, "docs", "archive_snapshot.tsv")
    open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("wrote %s for %s at %s" % (os.path.relpath(path, ROOT), tag, doi))


if __name__ == "__main__":
    main()
