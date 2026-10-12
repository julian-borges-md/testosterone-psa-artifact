"""Pass A: replay src/stage1_saturation/run1.py verbatim against the committed
union instrument set, then compare every cell against the committed
results/stage1/mr_results.tsv and mr_diagnostics.tsv.

Purpose: validate this environment end to end before any value is changed.
A cell that does not reproduce here is an environment problem, not a finding.
"""
import os, sys
import numpy as np
import pandas as pd

os.chdir("/tmp/gw")
sys.path.insert(0, "/tmp/gw")
import mr

REPO = sys.argv[1]
OUT = "/tmp/gw/out"
os.makedirs(OUT, exist_ok=True)

D = mr.load()
bim = pd.read_csv("union_snps.tsv", sep="\t", dtype={"chr": str}).set_index("rsid")
inst = {t: set(pd.read_csv(f"{t}_clump.clumped", sep=r"\s+").SNP)
        for t in ["biot", "shbg", "tt", "psa"]}


def region(s, chr_, lo, hi):
    r = bim.loc[s]
    return (r.chr == chr_) & (r.pos >= lo) & (r.pos <= hi)


def steiger(m):
    tx = m.b_x / m.se_x
    ty = m.b_y / m.se_y
    nx = m.n_x
    ny = m.n_y.fillna(m.n_y.median() if m.n_y.notna().any() else 300000)
    rx = tx ** 2 / (tx ** 2 + nx - 2)
    ry = ty ** 2 / (ty ** 2 + ny - 2)
    return rx > ry


R, M, H = [], [], []
plan = [("biot", "psa"), ("shbg", "psa"), ("tt", "psa"), ("biot", "prac"),
        ("biot", "fg"), ("shbg", "prac"), ("tt", "prac"), ("psa", "prac"),
        ("psa", "fg"), ("psa", "biot")]
for e, o in plan:
    ex = D[e][D[e].SNP.isin(inst[e])]
    lg = {}
    m = mr.harmonise(ex, D[o], lg)
    lg.update(exposure=e, outcome=o, instruments=len(ex))
    H.append(lg)
    lab = f"{e}->{o}"
    r, mt = mr.run(m, lab)
    R += r
    M.append(mt)
    st = steiger(m)
    r, mt = mr.run(m[st.values], lab + " | Steiger")
    R += r
    mt["steiger_removed"] = int((~st).sum())
    M.append(mt)
    excl = ~(region(m.SNP, "19", 50300000, 52400000).values
             | region(m.SNP, "17", 6500000, 8600000).values)
    if excl.sum() < len(m):
        r, mt = mr.run(m[excl], lab + " | excl SHBG/KLK3 loci")
        R += r
        mt["loci_removed"] = int((~excl).sum())
        M.append(mt)
    # keep the harmonised primary set for later reuse
    if (e, o) in [("biot", "psa"), ("tt", "psa"), ("shbg", "psa"), ("biot", "prac"),
                  ("tt", "prac"), ("shbg", "prac"), ("biot", "fg"), ("psa", "prac")]:
        m.to_csv(f"{OUT}/harm_{e}_{o}.tsv", sep="\t", index=False)
    print("done", lab, flush=True)

new = pd.DataFrame(R)
newm = pd.DataFrame(M)
new.to_csv(f"{OUT}/mr_results.tsv", sep="\t", index=False)
newm.to_csv(f"{OUT}/mr_diagnostics.tsv", sep="\t", index=False)
pd.DataFrame(H).to_csv(f"{OUT}/harmonisation_log.tsv", sep="\t", index=False)

# ---------------- reproduction check against committed ----------------
old = pd.read_csv(f"{REPO}/results/stage1/mr_results.tsv", sep="\t")
oldm = pd.read_csv(f"{REPO}/results/stage1/mr_diagnostics.tsv", sep="\t")
k = ["analysis", "method"]
j = old.merge(new, on=k, suffixes=("_old", "_new"))
j["d_beta"] = (j.beta_old - j.beta_new).abs()
j["d_se"] = (j.se_old - j.se_new).abs()
print("\n=== PASS A reproduction: %d committed cells, %d matched on key ==="
      % (len(old), len(j)))
print("max |d beta| = %.3e   max |d se| = %.3e" % (j.d_beta.max(), j.d_se.max()))
bad = j[(j.d_beta > 1e-9) | (j.d_se > 1e-9)]
print("cells differing beyond 1e-9: %d" % len(bad))
if len(bad):
    print(bad[["analysis", "method", "beta_old", "beta_new", "se_old", "se_new"]]
          .to_string(index=False))
jm = oldm.merge(newm, on="analysis", suffixes=("_old", "_new"))
for c in ["nsnp", "meanF", "Q", "egger_intercept"]:
    d = (jm[c + "_old"] - jm[c + "_new"]).abs().max()
    print("diagnostics %-16s max abs diff = %.3e" % (c, d))
print("\nharmonisation log (new):")
print(pd.DataFrame(H).to_string(index=False))
