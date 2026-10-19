import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
"""Write the regenerated estimator matrices into the repository as committed
result tables, one row per cell, so that every printed cell has a machine
readable source and the assertion module can check all of them.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

OUT = f"{GW}/out"
REPO = sys.argv[1]
DEST = os.path.join(REPO, "results", "stage7_matrix")
os.makedirs(DEST, exist_ok=True)

mrres = pd.read_csv(f"{OUT}/mr_results.tsv", sep="\t")
rob = {}
for f in glob.glob(f"{OUT}/robust_*.json"):
    r = json.load(open(f))
    rob[r["analysis"]] = r

ROWS = ["IVW", "MR Egger", "Weighted median", "Weighted mode",
        "IVW, Q outliers removed", "MR cML BIC (DP SE)", "MR RAPS Huber"]
ANALYSES = ["biot->psa", "tt->psa", "shbg->psa",
            "biot->prac", "tt->prac", "shbg->prac", "biot->fg", "psa->prac"]


def cell(analysis, row):
    if row == "IVW, Q outliers removed":
        d = mrres[(mrres.analysis == analysis)
                  & (mrres.method.str.startswith("IVW outliers removed"))].iloc[0]
        return d.beta, d.se, d.method
    if row == "MR cML BIC (DP SE)":
        return rob[analysis]["cml_dp_beta"], rob[analysis]["cml_dp_se"], \
            "MR cML BIC, data perturbation SE (B=500), K=%d" % rob[analysis]["cml_dp_K"]
    if row == "MR RAPS Huber":
        return rob[analysis]["raps_beta"], rob[analysis]["raps_se_boot"], \
            "MR RAPS Huber, bootstrap SE (B=300), tau2=%.6e" % rob[analysis]["raps_tau2"]
    d = mrres[(mrres.analysis == analysis) & (mrres.method == row)].iloc[0]
    return d.beta, d.se, d.method


rows = []
for a in ANALYSES:
    for r in ROWS:
        b, s, detail = cell(a, r)
        rows.append(dict(
            analysis=a, estimator=r, nsnp=rob[a]["nsnp"], beta=b, se=s,
            lo=b - 1.96 * s, hi=b + 1.96 * s, p=2 * stats.norm.sf(abs(b / s)),
            OR=np.exp(b), OR_lo=np.exp(b - 1.96 * s), OR_hi=np.exp(b + 1.96 * s),
            excludes_null=bool((b - 1.96 * s > 0) or (b + 1.96 * s < 0)),
            Q=rob[a]["Q"], I2=rob[a]["I2"], raps_tau2=rob[a]["raps_tau2"],
            detail=detail))
M = pd.DataFrame(rows)
M.to_csv(f"{DEST}/estimator_matrix.tsv", sep="\t", index=False)

# contrast table, all seven estimators
con = []
for r in ROWS:
    b1, s1, _ = cell("biot->psa", r)
    b2, s2, _ = cell("tt->psa", r)
    d = b1 - b2
    sd = np.sqrt(s1 ** 2 + s2 ** 2)
    con.append(dict(estimator=r, biot=b1, tt=b2, difference=d, se=sd, z=d / sd,
                    p=2 * stats.norm.sf(abs(d / sd)),
                    both_negative=bool(b1 < 0 and b2 < 0),
                    ratio=(b1 / b2 if abs(b2) > 1e-12 else np.nan)))
C = pd.DataFrame(con)
C.to_csv(f"{DEST}/contrast_biot_minus_tt.tsv", sep="\t", index=False)

# heterogeneity and overdispersion, per analysis
H = pd.DataFrame([dict(analysis=a, nsnp=rob[a]["nsnp"], Q=rob[a]["Q"], I2=rob[a]["I2"],
                       raps_tau2=rob[a]["raps_tau2"],
                       ivw_se_random=rob[a]["ivw_se"],
                       raps_se_bootstrap=rob[a]["raps_se_boot"],
                       cml_K=rob[a]["cml_K"], cml_dp_se=rob[a]["cml_dp_se"])
                  for a in ANALYSES])
H.to_csv(f"{DEST}/heterogeneity.tsv", sep="\t", index=False)

print("wrote %d matrix cells, %d contrast rows, %d heterogeneity rows to %s"
      % (len(M), len(C), len(H), DEST))
print(M.groupby("analysis").excludes_null.sum().to_string())
