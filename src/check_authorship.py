#!/usr/bin/env python3
"""Authorship and identifier integrity gate.

This repository is public, has a single author, and is the verification
chain for a manuscript. Three properties must hold. Two are closed by
construction, meaning that something nobody anticipated fails rather than
passing for want of being listed.

P1  Sole authorship in commit messages. Where the final paragraph of a
    message is a trailer block, only the keys in ALLOWED_TRAILER_KEYS may
    appear. Any other key fails, whatever name or address it carries.

P2  Closed URL world in commit messages. Only the hosts in ALLOWED_HOSTS
    may appear. A console, run or session link therefore cannot enter the
    public history by being a host nobody thought to forbid.

P3  No third party authorship attribution, product identification or
    authoring disclosure, in commit messages or in tracked documents. The
    forbidden terms are held as salted digests rather than plaintext, so
    this file does not itself contain the strings it exists to keep out of
    the repository. The technique is the one the archive citation gate
    uses for self referential matching.

    Class A covers attribution and product identification and is checked
    against every tracked text file and every commit message. Class B
    covers disclosure constructions and is checked against tracked
    documents only, so that a cited work on the subject is not mistaken
    for a disclosure of this work's authorship.

Findings can be excepted in docs/authorship_gate_exceptions.tsv, keyed by
commit sha or by repository relative path. Every exception is printed on
every run, so an accepted deviation stays visible instead of passing
silently.

Usage
  python src/check_authorship.py                 history and files
  python src/check_authorship.py --history-only [rev]
  python src/check_authorship.py --files-only
  python src/check_authorship.py --digest "a term"   print its digest
"""
import hashlib
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXCEPTIONS = ROOT / "docs" / "authorship_gate_exceptions.tsv"
SELF = pathlib.Path(__file__).resolve()

SALT = "tpsa-authorship-gate-v1"

# P1: the only trailer keys a commit message may carry.
ALLOWED_TRAILER_KEYS = {"signed-off-by"}

# P2: the only URL hosts a commit message may reference.
ALLOWED_HOSTS = {
    "doi.org", "dx.doi.org", "zenodo.org", "github.com",
    "www.ncbi.nlm.nih.gov", "pubmed.ncbi.nlm.nih.gov",
}

# P3: salted digests. See the module docstring for the two classes.
CLASS_A = {
    "0b06f829f59ec7fb",
    "10929784ff426409",
    "3a328ff7aba29537",
    "467a580ea4e0512a",
    "570e93672cd49ea7",
    "647a95b2f0ebc233",
    "6af1a008fb9519c1",
    "84bc8a9c1eca272a",
    "875fad8afa398356",
    "9274b9bb1f694d31",
    "99b1d77212f0de04",
    "a9d293f95cb67200",
    "aa543d8cd7b560b5",
    "b6185d72ebba4b7f",
    "bed2522f17f0a730",
    "c1a5c0e54baf555e",
    "cb94a2f9c7657901",
    "cbfce60091a81815",
    "d4c4ae5163c0fd69",
    "db63b31e80f9c6a2",
    "e178bdb16adb5109",
    "f718eb5f28f33a40"
}
CLASS_B = {
    "1559818d351ce83c",
    "4bf659ead7b1895a",
    "70fa0d85e3f75931",
    "80786710234f82ce",
    "98bba56733038d94",
    "9a8971add2cd476d",
    "9ec18be171554972",
    "a2ba09b78bf1ca22",
    "a5ab3f6923ef5cd5",
    "b4d35781b17f1923",
    "cb28126df96b3f21",
    "eb3c1e76b65a0047",
    "f5df812f9e0c2c44",
    "fef2f309fe33b9c4"
}
NGRAM_A = 5
NGRAM_B = 7

TEXT_SUFFIXES = {
    ".md", ".txt", ".py", ".r", ".sh", ".yml", ".yaml", ".cff",
    ".json", ".bib", ".tex", ".toml", ".cfg", ".ini",
}
TEXT_NAMES = {"Makefile", "LICENSE", "CITATION.cff", ".gitignore"}
DOCUMENT_SUFFIXES = {".md", ".cff", ".txt"}
SKIP = {"docs/authorship_gate_exceptions.tsv", "src/check_authorship.py"}
MAX_BYTES = 2_000_000

TRAILER_LINE = re.compile(r"^([A-Za-z][A-Za-z0-9_-]*):\s")
URL = re.compile(r"https?://([A-Za-z0-9._~-]+)", re.IGNORECASE)
WORDS = re.compile(r"[a-z0-9]+")


def digest(term):
    norm = " ".join(term.lower().split())
    return hashlib.sha256((SALT + "|" + norm).encode()).hexdigest()[:16]


def term_hits(text, digests, span):
    """Return the digests present in text, by hashing word n-grams."""
    if not digests:
        return set()
    toks = WORDS.findall(text.lower())
    found = set()
    for n in range(1, span + 1):
        if n > len(toks):
            break
        for i in range(len(toks) - n + 1):
            h = digest(" ".join(toks[i:i + n]))
            if h in digests:
                found.add(h)
    return found


def check_message(msg):
    hits = []
    paras = [p for p in msg.strip().split("\n\n") if p.strip()]
    if paras:
        last = [l for l in paras[-1].splitlines() if l.strip()]
        keyed = [TRAILER_LINE.match(l) for l in last]
        if last and all(keyed):
            for m in keyed:
                key = m.group(1).lower()
                if key not in ALLOWED_TRAILER_KEYS:
                    hits.append("P1 trailer key not permitted: " + m.group(1))
    for m in URL.finditer(msg):
        host = m.group(1).lower()
        if host not in ALLOWED_HOSTS:
            hits.append("P2 URL host not permitted: " + host)
    for h in sorted(term_hits(msg, CLASS_A, NGRAM_A)):
        hits.append("P3 class A term present, digest " + h)
    return hits


def check_file(rel, text):
    hits = []
    for h in sorted(term_hits(text, CLASS_A, NGRAM_A)):
        hits.append("P3 class A term present, digest " + h)
    if pathlib.Path(rel).suffix.lower() in DOCUMENT_SUFFIXES:
        for h in sorted(term_hits(text, CLASS_B, NGRAM_B)):
            hits.append("P3 class B term present, digest " + h)
    return hits


def read_exceptions():
    out = {}
    if not EXCEPTIONS.exists():
        return out
    for line in EXCEPTIONS.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        if parts and parts[0] and parts[0] != "key":
            out[parts[0]] = " ".join(parts[1:])
    return out


def git(*args):
    return subprocess.run(["git", "-C", str(ROOT), *args],
                          capture_output=True, text=True, check=True).stdout


def is_text(rel):
    p = pathlib.Path(rel)
    return p.suffix.lower() in TEXT_SUFFIXES or p.name in TEXT_NAMES


def main(argv):
    if "--digest" in argv:
        i = argv.index("--digest")
        print(digest(argv[i + 1]))
        return 0

    do_history = "--files-only" not in argv
    do_files = "--history-only" not in argv
    rev = "HEAD"
    for a in argv[1:]:
        if not a.startswith("--") and argv[argv.index(a) - 1] != "--digest":
            rev = a

    exceptions = read_exceptions()
    failures, excused = [], []
    n_commits = n_files = 0

    if do_history:
        for sha in git("rev-list", rev).split():
            n_commits += 1
            hits = check_message(git("log", "-1", "--format=%B", sha))
            if not hits:
                continue
            key = next((k for k in (sha, sha[:7]) if k in exceptions), None)
            (excused if key else failures).append(
                (sha[:7], hits, exceptions.get(key, "")))

    if do_files:
        for rel in git("ls-files").split("\n"):
            if not rel or rel in SKIP or not is_text(rel):
                continue
            path = ROOT / rel
            if not path.exists() or path.resolve() == SELF:
                continue
            if path.stat().st_size > MAX_BYTES:
                continue
            n_files += 1
            hits = check_file(rel, path.read_text(encoding="utf-8",
                                                 errors="ignore"))
            if not hits:
                continue
            (excused if rel in exceptions else failures).append(
                (rel, hits, exceptions.get(rel, "")))

    print("authorship gate: %d commits and %d tracked files checked"
          % (n_commits, n_files))
    for subj, hits, reason in excused:
        print("EXCEPTED  %s  %s" % (subj, reason))
        for h in hits:
            print("          " + h)
    if not failures:
        print("PASS  %d excepted, 0 failed" % len(excused))
        return 0
    for subj, hits, _ in failures:
        print("FAIL  " + subj)
        for h in hits:
            print("      " + h)
    print("FAIL  %d subjects violate the gate" % len(failures))
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
