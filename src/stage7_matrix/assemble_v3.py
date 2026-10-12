"""Assemble every number reported in the manuscript from the regenerated run, and compute
the derived quantities the R2 review requires. Prints a single auditable block.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammaln, logsumexp

OUT = "/tmp/gw/out"
REPO = sys.argv[1]
pd.set_option("display.width", 250)

mrres = pd.read_csv(f"{OUT}/mr_results.tsv", sep="\t")
diag = pd.read_csv(f"{OUT}/mr_diagnostics.tsv", sep="\t")
rob = {}
for f in glob.glob(f"{OUT}/robust_*.json"):
    r = json.load(open(f))
    rob[r["analysis"]] = r


def g(analysis, method):
    d = mrres[(mrres.analysis == analysis) & (mrres.method == method)]
    if len(d) == 0:
        d = mrres[(mrres.analysis == analysis)
                  & (mrres.method.str.startswith("IVW outliers removed"))]
    return d.iloc[0]


ROWS = ["IVW", "MR Egger", "Weighted median", "Weighted mode",
        "IVW, Q outliers removed", "MR cML BIC (DP SE)", "MR RAPS Huber"]


def cell(analysis, row):
    if row == "IVW, Q outliers removed":
        r = g(analysis, "QOUT")
        return r.beta, r.se
    if row == "MR cML BIC (DP SE)":
        return rob[analysis]["cml_dp_beta"], rob[analysis]["cml_dp_se"]
    if row == "MR RAPS Huber":
        return rob[analysis]["raps_beta"], rob[analysis]["raps_se_boot"]
    r = g(analysis, row)
    return r.beta, r.se


def fmt(b, s, scale):
    lo, hi = b - 1.96 * s, b + 1.96 * s
    sig = (lo > 0) or (hi < 0)
    if scale == "log":
        return "%.3f (%.3f to %.3f)%s" % (b, lo, hi, "  *" if sig else "")
    return "%.3f (%.3f to %.3f)%s" % (np.exp(b), np.exp(lo), np.exp(hi), "  *" if sig else "")


print("=" * 110)
print("TABLE 2  effect on PSA, natural log units per SD of exposure   (* = 95% interval excludes the null)")
print("=" * 110)
cols2 = [("biot->psa", "Bioavailable T"), ("tt->psa", "Total T"), ("shbg->psa", "SHBG")]
t2 = []
for row in ROWS:
    rec = {"Estimator": row}
    for a, lab in cols2:
        b, s = cell(a, row)
        rec[lab] = fmt(b, s, "log")
    t2.append(rec)
rec = {"Estimator": "Instruments analysed"}
for a, lab in cols2:
    rec[lab] = str(rob[a]["nsnp"])
t2.append(rec)
rec = {"Estimator": "Cochran Q (I squared)"}
for a, lab in cols2:
    rec[lab] = "%.0f (%.1f%%)" % (rob[a]["Q"], rob[a]["I2"])
t2.append(rec)
rec = {"Estimator": "MR RAPS tau squared"}
for a, lab in cols2:
    rec[lab] = "%.2e" % rob[a]["raps_tau2"]
t2.append(rec)
print(pd.DataFrame(t2).to_string(index=False))

print()
print("=" * 110)
print("TABLE 4  prostate cancer, odds ratio per SD of exposure")
print("=" * 110)
cols4 = [("biot->prac", "Bioav T, PRACTICAL"), ("tt->prac", "Total T, PRACTICAL"),
         ("shbg->prac", "SHBG, PRACTICAL"), ("biot->fg", "Bioav T, FinnGen")]
t4 = []
for row in ROWS:
    rec = {"Estimator": row}
    for a, lab in cols4:
        b, s = cell(a, row)
        rec[lab] = fmt(b, s, "or")
    t4.append(rec)
rec = {"Estimator": "Instruments analysed"}
for a, lab in cols4:
    rec[lab] = str(rob[a]["nsnp"])
t4.append(rec)
rec = {"Estimator": "Cochran Q (I squared)"}
for a, lab in cols4:
    rec[lab] = "%.0f (%.1f%%)" % (rob[a]["Q"], rob[a]["I2"])
t4.append(rec)
rec = {"Estimator": "MR RAPS tau squared"}
for a, lab in cols4:
    rec[lab] = "%.2e" % rob[a]["raps_tau2"]
t4.append(rec)
print(pd.DataFrame(t4).to_string(index=False))

print()
print("=" * 110)
print("TABLE 3  contrast bioavailable minus total testosterone on PSA, ALL SEVEN estimators")
print("=" * 110)
t3 = []
for row in ROWS:
    b1, s1 = cell("biot->psa", row)
    b2, s2 = cell("tt->psa", row)
    d = b1 - b2
    sd = np.sqrt(s1 ** 2 + s2 ** 2)
    z = d / sd
    p = 2 * stats.norm.sf(abs(z))
    ratio = b1 / b2 if abs(b2) > 1e-12 else np.nan
    t3.append(dict(Estimator=row, bioav="%.4f" % b1, total="%.4f" % b2,
                   difference="%+.4f" % d, SE="%.4f" % sd, z="%+.2f" % z,
                   p=("%.4f" % p if p >= 1e-4 else "<0.0001"),
                   ratio="%.2f" % ratio,
                   note="both negative" if (b1 < 0 and b2 < 0) else ""))
T3 = pd.DataFrame(t3)
print(T3.to_string(index=False))
pos = [float(r["ratio"]) for _, r in T3.iterrows() if r["note"] == ""]
print("\nM2 ratio range across estimators where both estimates are positive: %.2f to %.2f fold"
      % (min(pos), max(pos)))
print("M2 note: the IVW denominator is %.4f (SE %.4f), indistinguishable from zero, so the ratio is unstable."
      % cell("tt->psa", "IVW"))

print()
print("=" * 110)
print("DECOMPOSITION  (M5)")
print("=" * 110)
mv = pd.read_csv(f"{REPO}/results/stage1/mvmr_results.tsv", sep="\t")
d_prac = mv[(mv.model == "biot+psa->prac") & (mv.exposure == "biot")].iloc[0]
p_prac = mv[(mv.model == "biot+psa->prac") & (mv.exposure == "psa")].iloc[0]
d_fg = mv[(mv.model == "biot+psa->fg") & (mv.exposure == "biot")].iloc[0]
a = g("biot->psa", "IVW").beta
tot = g("biot->prac", "IVW").beta
indirect = a * p_prac.beta
prod_share = indirect / tot
diff_share = (tot - d_prac.beta) / tot
gap = tot - (indirect + d_prac.beta)
print("total effect (IVW)                  %.6f" % tot)
print("direct effect (MVMR, PSA held)      %.6f   OR %.3f (%.3f to %.3f)"
      % (d_prac.beta, np.exp(d_prac.beta), np.exp(d_prac.lo), np.exp(d_prac.hi)))
print("PSA on cancer, testosterone held    %.6f" % p_prac.beta)
print("indirect, product of coefficients   %.6f" % indirect)
print("channel share, product method       %.4f  (%.1f%%)" % (prod_share, 100 * prod_share))
print("channel share, difference method    %.4f  (%.1f%%)" % (diff_share, 100 * diff_share))
print("ratio of the two methods            %.3f" % (diff_share / prod_share))
print("additivity gap  total-(ind+dir)     %.6f  (%.1f%% of the total effect)"
      % (gap, 100 * gap / tot))
print("conditional F, biot+psa->PRACTICAL  %.2f" % d_prac.condF)
print("conditional F, biot+psa->FinnGen    %.2f" % d_fg.condF)
print("case fraction PRACTICAL             %.3f  (79148 cases, 61106 controls)"
      % (79148 / (79148 + 61106)))

st_psa = g("psa->prac | Steiger", "IVW").beta
un_psa = g("psa->prac", "IVW").beta
st_rm = int(diag[diag.analysis == "psa->prac | Steiger"].iloc[0].steiger_removed)
n_psa = int(diag[diag.analysis == "psa->prac"].iloc[0].nsnp)
print("\nSteiger filtering, PSA on cancer    %d of %d instruments removed" % (st_rm, n_psa))
print("PSA on cancer, unfiltered / filtered %.4f / %.4f  (%.1f%% reduction)"
      % (un_psa, st_psa, 100 * (1 - st_psa / un_psa)))
print("channel share under Steiger filter   %.4f  (%.1f%%)"
      % (a * st_psa / tot, 100 * a * st_psa / tot))
print("channel share range across specifications  %.1f%% to %.1f%%"
      % (100 * a * st_psa / tot, 100 * prod_share))

print()
print("=" * 110)
print("SATURATION BOUND  (M4, Mod1)")
print("=" * 110)
# Vermeulen 1999, calibrated against the published UK Biobank male free testosterone mean
Kt, Ka, MWalb = 1.0e9, 3.6e4, 66500.0


def vermeulen(T_nM, S_nM, alb_gL):
    T, S = T_nM * 1e-9, S_nM * 1e-9
    A = alb_gL / MWalb
    N = 1 + Ka * A
    aa = N * Kt
    bb = N + Kt * (S - T)
    FT = (-bb + np.sqrt(bb ** 2 + 4 * aa * T)) / (2 * aa)
    return FT * 1e12, N * FT * 1e9      # free in pmol/L, bioavailable in nmol/L


TT_MEAN, TT_SD = 12.2, 3.5               # Hynes 2024, UK Biobank men, n = 16,237
SHBG_MED, ALB = 38.9, 43.0
ft, bt = vermeulen(TT_MEAN, SHBG_MED, ALB)
print("calibration: total T %.1f, SHBG %.1f, albumin %.0f  ->  free T %.1f pmol/L (published 214.4), "
      "bioavailable T %.2f nmol/L" % (TT_MEAN, SHBG_MED, ALB, ft, bt))

sig = np.sqrt(np.log(1 + (TT_SD / TT_MEAN) ** 2))
mu = np.log(TT_MEAN) - 0.5 * sig ** 2
print("lognormal fit to total testosterone: mu %.5f  sigma %.5f" % (mu, sig))
print("\n thresh_totalT  bioavail_equiv   frac_below   required_slope_ratio (1/frac)")
for thr in [4.0, 6.0, 8.0]:
    _, bthr = vermeulen(thr, SHBG_MED, ALB)
    frac = stats.norm.cdf((np.log(thr) - mu) / sig)
    print("   %4.1f nmol/L    %5.2f nmol/L   %8.4f%%   %12.1f" % (thr, bthr, 100 * frac, 1 / frac))
print("\nPreviously used distribution (mean 11.9, SD 3.6), for comparison:")
sig0 = np.sqrt(np.log(1 + (3.6 / 11.9) ** 2))
mu0 = np.log(11.9) - 0.5 * sig0 ** 2
for thr in [4.0, 6.0, 8.0]:
    frac = stats.norm.cdf((np.log(thr) - mu0) / sig0)
    print("   %4.1f nmol/L                     %8.4f%%   %12.1f" % (thr, 100 * frac, 1 / frac))

# sensitivity of the bioavailable threshold to SHBG
print("\nbioavailable equivalent of the 8 nmol/L total testosterone threshold, across SHBG:")
for s in [20, 30, 38.9, 50, 70]:
    _, b8 = vermeulen(8.0, s, ALB)
    print("   SHBG %4.1f nmol/L  ->  %.2f nmol/L bioavailable" % (s, b8))

print()
print("=" * 110)
print("COMPATIBILITY WITH RANDOMIZED EVIDENCE  (M6)")
print("=" * 110)
P = pd.read_csv(f"{REPO}/results/stage4/pooled_priors.tsv", sep="\t")
tot_p = P[P.prior == "total pooled"].iloc[0]
se_zeng = (np.log(1.62) - np.log(0.44)) / 3.92
for lab, b, s in [("pooled prior vs Zeng RCT only", tot_p.beta, tot_p.se),
                  ("PRACTICAL only prior vs Zeng RCT only", 0.1762, 0.0511)]:
    z = (b - np.log(0.85)) / np.sqrt(s ** 2 + se_zeng ** 2)
    print("%-42s z %.2f  p %.3f" % (lab, z, 2 * stats.norm.sf(abs(z))))
print("committed results/stage3/conflict.txt, 41 trial likelihood: 1 SD p 0.272, 2 SD p 0.098")
print("committed results/stage3/zeng_rct_only.txt, RR 0.85 anchor: 1 SD p 0.314, 2 SD p 0.139")
print("committed zeng_rct_only.txt, RR 0.62 anchor: 1 SD p 0.007, 2 SD p 0.001  "
      "<- anchor WITHDRAWN: pools non randomized cohorts")
print("live face value compatibility range: p 0.098 to 0.369")

print()
print("=" * 110)
print("TRIAL LIKELIHOOD AND REVERSE BAYES  (Mod3, M6)")
print("=" * 110)
T = pd.read_csv(f"{REPO}/results/stage3/trial_counts_from_sharifan_supp.tsv", sep="\t")


def lc(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


gr = np.linspace(np.log(0.05), np.log(20), 4001)


def loglik(df):
    ll = np.zeros_like(gr)
    for _, t in df.iterrows():
        m = t.e1 + t.e0
        xs = np.arange(max(0, m - t.n0), min(m, t.n1) + 1)
        ll += lc(t.n1, t.e1) + lc(t.n0, t.e0) + t.e1 * gr - logsumexp(
            (lc(t.n1, xs) + lc(t.n0, m - xs))[None, :] + np.outer(gr, xs), axis=1)
    return ll - ll.max()


ll = loglik(T)
ci = gr[ll >= -1.92]
sd_lik = (ci.max() - ci.min()) / 3.92
print("trials entering the likelihood: %d, participants %d, events %d (largest trial supplies %d)"
      % (len(T), T.n1.sum() + T.n0.sum(), T.e1.sum() + T.e0.sum(), (T.e1 + T.e0).max()))
print("exact conditional OR %.3f (profile 95%% CI %.3f to %.3f)"
      % (np.exp(gr[ll.argmax()]), np.exp(ci.min()), np.exp(ci.max())))
print("trial likelihood SD on the log odds scale: %.4f   <- the firmness denominator" % sd_lik)


def post(mu0, s0, llv=ll):
    lp = llv + stats.norm.logpdf(gr, mu0, s0)
    p = np.exp(lp - lp.max())
    p /= p.sum()
    c = np.cumsum(p)
    return (np.exp(gr[np.searchsorted(c, 0.5)]), np.exp(gr[np.searchsorted(c, 0.025)]),
            np.exp(gr[np.searchsorted(c, 0.975)]), p[gr > 0].sum(), p[gr > np.log(1.2)].sum())


for k, lab in [(0.5, "twice as firm (prior SD = half the likelihood SD)"),
               (1.0, "as firm as the trials (prior SD = likelihood SD)"),
               (2.0, "half as firm (prior SD = twice the likelihood SD)")]:
    s0 = sd_lik * k
    lo, hi = -2.0, 5.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if post(mid, s0)[3] >= 0.95:
            hi = mid
        else:
            lo = mid
    print("  prior SD %.3f  %-52s minimum prior median OR %.2f" % (s0, lab, np.exp(hi)))

# leave one out on the largest trial
T2 = T[(T.e1 + T.e0) < (T.e1 + T.e0).max()]
ll2 = loglik(T2)
dr = P[P.prior == "direct pooled"].iloc[0]
print("\nleave one out, excluding the largest trial (%s, %d events):"
      % (T.loc[(T.e1 + T.e0).idxmax(), "trial"], (T.e1 + T.e0).max()))
print("  exact conditional OR %.3f ; PSA removed posterior OR %.3f (%.3f to %.3f), P(OR>1.2) %.3f"
      % (np.exp(gr[ll2.argmax()]), *post(dr.beta, dr.se, ll2)[:3], post(dr.beta, dr.se, ll2)[4]))
print("  full set,            PSA removed posterior OR %.3f (%.3f to %.3f), P(OR>1.2) %.3f"
      % (*post(dr.beta, dr.se)[:3], post(dr.beta, dr.se)[4]))

print("\nposterior across assumed therapeutic shifts, PSA channel removed, pooled prior:")
vals = []
for k in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    m = post(k * dr.beta, k * dr.se)
    vals.append(m[0])
    print("   shift %.1f SD  posterior OR %.4f" % (k, m[0]))
print("   range across 0.5 to 3 SD: %.3f to %.3f"
      % (min(vals), max(vals)))

print()
print("=" * 110)
print("POWER TABLE ASSUMPTIONS  (Mod3)")
print("=" * 110)
pw = pd.read_csv(f"{REPO}/results/stage2/power.tsv", sep="\t")
sc = pw[pw.scenario == "1.19"][["cases", "case_frac"]].copy()
sc["controls"] = (sc.cases / sc.case_frac - sc.cases).round(0).astype(int)
sc["total"] = (sc.cases / sc.case_frac).round(0).astype(int)
print(sc.to_string(index=False))
print("controls are held between %d and %d across every scenario; R squared 0.040; alpha 0.05"
      % (sc.controls.min(), sc.controls.max()))
tost = pw[pw.scenario == "TOST(true OR=1)"][["cases", "value"]]
print("two one sided test power for equivalence within OR 1.10:")
print(tost.to_string(index=False))
