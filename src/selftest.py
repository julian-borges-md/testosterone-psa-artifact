"""Proves the assertion modules can fail.

Two perturbations, each in a temporary copy so the tree is never touched.
A module that cannot fail asserts nothing, so every module whose output is
quoted as evidence needs a case here.

  verify.py           one committed estimate is altered
  check_manuscript.py one count on the title page is altered
"""
import io, os, re, shutil, subprocess, sys, tempfile
import pandas as pd

FAILURES = []


def expect_nonzero(label, argv, cwd=None, saying=None):
    """Require a non-zero exit, and that it is for the perturbation's own reason.

    Exit status alone is not evidence: a module with any other outstanding
    discordance exits non-zero whatever is perturbed, which would make the case
    vacuous. When `saying` is given, the perturbation's own message must appear.
    """
    r = subprocess.run(argv, capture_output=True, text=True, cwd=cwd)
    out = (r.stdout or "") + (r.stderr or "")
    if r.returncode == 0:
        FAILURES.append("%s: the module did not detect the perturbation" % label)
    elif saying and saying not in out:
        FAILURES.append("%s: exited %d but never reported %r, so the exit is "
                        "not evidence that the perturbation was detected"
                        % (label, r.returncode, saying))
    else:
        print("Selftest passed: %s detected the perturbation and exited %d"
              % (label, r.returncode))


# -- verify.py: a committed estimate is altered -----------------------------
tmp = tempfile.mkdtemp()
try:
    shutil.copytree("results", os.path.join(tmp, "results"))
    p = os.path.join(tmp, "results/stage1/mr_results.tsv")
    d = pd.read_csv(p, sep="\t")
    d.loc[(d.analysis == "biot->psa") & (d.method == "IVW"), "beta"] = 0.0
    d.to_csv(p, sep="\t", index=False)
    expect_nonzero("verify.py",
                   [sys.executable, "src/verify.py", os.path.join(tmp, "results/")])
finally:
    shutil.rmtree(tmp, ignore_errors=True)

# -- check_manuscript.py: a title page count is altered ---------------------
# The title page claims its counts cannot drift from the manuscript. That claim
# is only true while this case passes.
tmp = tempfile.mkdtemp()
try:
    for d in ("manuscript", "results", "docs", "figures", "src"):
        if os.path.isdir(d):
            shutil.copytree(d, os.path.join(tmp, d))
    for f in ("CITATION.cff", "README.md"):
        if os.path.isfile(f):
            shutil.copy(f, os.path.join(tmp, f))
    tp = os.path.join(tmp, "manuscript", "title_page.md")
    body = io.open(tp, encoding="utf-8").read()
    moved = re.sub(r"(Main text )([\d,]+)( words)",
                   lambda m: m.group(1) + "9,999" + m.group(3), body, count=1)
    if moved == body:
        FAILURES.append("check_manuscript.py: the title page states no main "
                        "text count, so there was nothing to perturb")
    else:
        io.open(tp, "w", encoding="utf-8").write(moved)
        expect_nonzero("check_manuscript.py",
                       [sys.executable, "src/check_manuscript.py",
                        "manuscript/manuscript_v6.0.md"], cwd=tmp,
                       saying="title page: states 9,999 for main text words")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

# -- check_manuscript.py: an unbounded share claim --------------------------
# The gate is inert on a correct manuscript, so a pass proves nothing about it.
# This case restores the wording the gate exists to catch.
tmp = tempfile.mkdtemp()
try:
    for d in ("manuscript", "results", "docs", "figures", "src"):
        if os.path.isdir(d):
            shutil.copytree(d, os.path.join(tmp, d))
    for f in ("CITATION.cff", "README.md"):
        if os.path.isfile(f):
            shutil.copy(f, os.path.join(tmp, f))
    ms = os.path.join(tmp, "manuscript", "manuscript_v6.0.md")
    body = io.open(ms, encoding="utf-8").read()
    moved = body.replace("part of its apparent cancer risk",
                         "much of its apparent cancer risk", 1)
    if moved == body:
        FAILURES.append("check_manuscript.py: the conclusion no longer carries "
                        "the bounded share wording, so there was nothing to "
                        "perturb; update this case to the current wording")
    else:
        io.open(ms, "w", encoding="utf-8").write(moved)
        expect_nonzero("check_manuscript.py share gate",
                       [sys.executable, "src/check_manuscript.py",
                        "manuscript/manuscript_v6.0.md"], cwd=tmp,
                       saying="share claim: a magnitude is asserted")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if FAILURES:
    print("\nSELFTEST FAILED:", file=sys.stderr)
    for f in FAILURES:
        print("  " + f, file=sys.stderr)
    sys.exit(1)
