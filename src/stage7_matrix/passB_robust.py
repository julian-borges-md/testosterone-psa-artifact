"""Pass B: regenerate the robust estimator rows for every cell of the estimator
matrix from one run.

Two properties this module guarantees:

  1. MR RAPS. The committed estimator (src/stage2_detection/robust2.py, raps())
     returns a 300 replicate nonparametric bootstrap standard error and the
     fitted overdispersion parameter tau squared. A non overdispersed fit reports standard errors
     equal to the fixed effect IVW values, which removes overdispersion. This
     script restores the committed estimator and reports tau squared alongside
     every estimate.

  2. Coverage. The committed run covered six exposure to outcome pairs, so four
     cells of Table 4 (MR cML and MR RAPS for total testosterone and
     for SHBG on PRACTICAL) had no committed provenance at all. This script
     covers all eight pairs needed by Tables 2, 3 and 4.

Change from the committed run: the random number
generator is reseeded per analysis, so each cell is independent of the order in
which the others are computed. The first analysis in the committed order
(biot->psa) therefore reproduces bit for bit; later cells differ from the
committed values only by bootstrap Monte Carlo noise.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import optimize, stats

os.chdir("/tmp/gw")
sys.path.insert(0, "/tmp/gw")
import mr

OUT = "/tmp/gw/out"
os.makedirs(OUT, exist_ok=True)
SEED_MAIN = 20261011
SEED_DP = 7

rng = np.random.default_rng(SEED_MAIN)


# ---------- MR cML (Xue, Shen, Pan 2021): committed implementation ----------
def cml_fit(bx, by, sx, sy, K, theta0, iters=200):
    th = theta0
    n = len(bx)
    r = np.zeros(n)
    for _ in range(iters):
        bj = (bx / sx ** 2 + th * (by - r) / sy ** 2) / (1 / sx ** 2 + th ** 2 / sy ** 2)
        res = (by - th * bj) ** 2 / sy ** 2
        inv = np.argsort(-res)[:K] if K > 0 else np.array([], int)
        r = np.zeros(n)
        r[inv] = by[inv] - th * bj[inv]
        v = np.ones(n, bool)
        v[inv] = False
        th_new = np.sum(bj[v] * by[v] / sy[v] ** 2) / np.sum(bj[v] ** 2 / sy[v] ** 2)
        if abs(th_new - th) < 1e-10:
            th = th_new
            break
        th = th_new
    bj = (bx / sx ** 2 + th * (by - r) / sy ** 2) / (1 / sx ** 2 + th ** 2 / sy ** 2)
    nll = 0.5 * np.sum((bx - bj) ** 2 / sx ** 2) + 0.5 * np.sum((by - th * bj - r) ** 2 / sy ** 2)
    v = r == 0
    se = 1 / np.sqrt(np.sum(bx[v] ** 2 / (sy[v] ** 2 + th ** 2 * sx[v] ** 2)))
    return th, se, nll, v


def cml(bx, by, sx, sy, nmin, starts=10):
    n = len(bx)
    fits = []
    init = np.sum(bx * by / sy ** 2) / np.sum(bx ** 2 / sy ** 2)
    for K in range(0, n - 2):
        best = None
        for s in np.concatenate([[init], rng.uniform(init - 0.5, init + 0.5, starts)]):
            f = cml_fit(bx, by, sx, sy, K, s)
            if best is None or f[2] < best[2]:
                best = f
        bic = 2 * best[2] + np.log(nmin) * K
        fits.append((K, best[0], best[1], bic, best[3]))
        if K > n // 2:
            break
    F = pd.DataFrame(fits, columns=["K", "theta", "se", "bic", "valid"])
    b = F.loc[F.bic.idxmin()]
    w = np.exp(-(F.bic - F.bic.min()) / 2)
    w /= w.sum()
    ma = np.sum(w * F.theta)
    mase = np.sqrt(np.sum(w * np.sqrt(F.se ** 2 + (F.theta - ma) ** 2)) ** 2)
    return dict(K_bic=int(b.K), theta_bic=b.theta, se_bic=b.se, theta_ma=ma, se_ma=mase)


# ---------- MR RAPS with Huber loss (Zhao 2020): committed implementation ----------
def raps_est(bx, by, sx, sy, c=1.345):
    def obj_full(p):
        th, lt = p
        tau2 = np.exp(lt)
        V = sy ** 2 + th ** 2 * sx ** 2 + tau2
        t = (by - th * bx) / np.sqrt(V)
        a = np.abs(t)
        rho = np.where(a <= c, 0.5 * t ** 2, c * a - 0.5 * c ** 2)
        return np.sum(rho) + 0.5 * np.sum(np.log(V))

    init = np.sum(bx * by / sy ** 2) / np.sum(bx ** 2 / sy ** 2)
    r = optimize.minimize(obj_full, [init, np.log(np.median(sy ** 2))],
                          method="Nelder-Mead",
                          options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
    return r.x[0], np.exp(r.x[1])


def raps(bx, by, sx, sy, B=300):
    th, tau2 = raps_est(bx, by, sx, sy)
    bs = []
    for _ in range(B):
        i = rng.integers(0, len(bx), len(bx))
        bs.append(raps_est(bx[i], by[i], sx[i], sy[i])[0])
    return th, np.std(bs), tau2


# ---------- cML data perturbation standard error ----------
def cml_fast(bx, by, sx, sy, nmin, Kmax):
    best = None
    init = np.sum(bx * by / sy ** 2) / np.sum(bx ** 2 / sy ** 2)
    for K in range(0, Kmax + 1):
        f = min((cml_fit(bx, by, sx, sy, K, s, iters=100) for s in [init, init - 0.2, init + 0.2]),
                key=lambda z: z[2])
        bic = 2 * f[2] + np.log(nmin) * K
        if best is None or bic < best[1]:
            best = (f[0], bic, K)
    return best[0], best[2]


PLAN = [
    ("biot", "psa", 178782),
    ("tt", "psa", 194453),
    ("shbg", "psa", 180726),
    ("biot", "prac", 140254),
    ("tt", "prac", 140254),
    ("shbg", "prac", 140254),
    ("biot", "fg", 178782),
    ("psa", "prac", 140254),
]

only = sys.argv[1:] if len(sys.argv) > 1 else None

for e, o, nmin in PLAN:
    tag = f"{e}_{o}"
    if only and tag not in only:
        continue
    dest = f"{OUT}/robust_{tag}.json"
    if os.path.exists(dest):
        print("skip (exists)", tag, flush=True)
        continue
    m = pd.read_csv(f"{OUT}/harm_{tag}.tsv", sep="\t")
    bx, by, sx, sy = m.b_x.values, m.b_y.values, m.se_x.values, m.se_y.values
    n = len(m)

    # reseed so this cell does not depend on the order of the others
    rng = np.random.default_rng(SEED_MAIN)
    globals()["rng"] = rng

    ivw_b, ivw_s, Q, k = mr.ivw(bx, by, sy)
    i2 = max(0.0, (Q - (k - 1)) / Q) * 100

    c = cml(bx, by, sx, sy, nmin)
    rt, rse, tau2 = raps(bx, by, sx, sy)

    # data perturbation standard error for the cML BIC point estimate
    Kmax = int(min(n // 2, max(70, c["K_bic"] + 20)))
    dprng = np.random.default_rng(SEED_DP)
    est, Ksel = cml_fast(bx, by, sx, sy, nmin, Kmax)
    bs = [cml_fast(bx + dprng.normal(0, sx), by + dprng.normal(0, sy), sx, sy, nmin, Kmax)[0]
          for _ in range(500)]
    dp_se = float(np.std(bs))

    rec = dict(
        analysis=f"{e}->{o}", nsnp=int(n), Q=float(Q), I2=float(i2),
        ivw_beta=float(ivw_b), ivw_se=float(ivw_s),
        cml_K=int(c["K_bic"]), cml_beta=float(c["theta_bic"]), cml_se_profile=float(c["se_bic"]),
        cml_ma_beta=float(c["theta_ma"]), cml_ma_se=float(c["se_ma"]),
        cml_dp_beta=float(est), cml_dp_K=int(Ksel), cml_dp_Kmax=Kmax, cml_dp_se=dp_se,
        raps_beta=float(rt), raps_se_boot=float(rse), raps_tau2=float(tau2),
        raps_B=300, dp_B=500,
    )
    with open(dest, "w") as fh:
        json.dump(rec, fh, indent=1)
    print(json.dumps(rec), flush=True)

print("pass B complete", flush=True)
