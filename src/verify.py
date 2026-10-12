"""RLO_TPC | Testosterone, PSA and prostate cancer
Re-derives every number asserted in the manuscript from the committed result
tables and the committed trial counts. No network access and no raw GWAS
download required. Exits non-zero and prints DRIFT to stderr on any mismatch.
"""
import sys, numpy as np, pandas as pd
from scipy import stats, optimize
from scipy.special import gammaln, logsumexp

R = sys.argv[1] if len(sys.argv) > 1 else 'results/'
fails, rows = [], []

def check(cid, label, got, expected, tol):
    ok = abs(got - expected) <= tol
    rows.append((cid, label, got, expected, 'ok' if ok else 'DRIFT'))
    if not ok: fails.append(cid)

def mr(path, analysis, method):
    d = pd.read_csv(R + path, sep='\t')
    return d[(d.analysis == analysis) & (d.method == method)].iloc[0]

# ---------------- Stage 1: saturation bioassay ----------------
check('F1a', 'Bioavailable T on PSA, IVW', mr('stage1/mr_results.tsv', 'biot->psa', 'IVW').beta, 0.085, 6e-4)
check('F1b', 'Bioavailable T on PSA, weighted median', mr('stage1/mr_results.tsv', 'biot->psa', 'Weighted median').beta, 0.130, 6e-4)
check('F1c', 'Bioavailable T on PSA, weighted mode', mr('stage1/mr_results.tsv', 'biot->psa', 'Weighted mode').beta, 0.149, 6e-4)
check('F2a', 'Total T on PSA, IVW', mr('stage1/mr_results.tsv', 'tt->psa', 'IVW').beta, 0.010, 6e-4)
check('F2b', 'SHBG on PSA, IVW', mr('stage1/mr_results.tsv', 'shbg->psa', 'IVW').beta, -0.076, 6e-4)
check('F3a', 'Bioavailable T on prostate cancer, PRACTICAL, OR', np.exp(mr('stage1/mr_results.tsv', 'biot->prac', 'IVW').beta), 1.19, 6e-3)
check('F3b', 'Bioavailable T on prostate cancer, FinnGen, OR', np.exp(mr('stage1/mr_results.tsv', 'biot->fg', 'IVW').beta), 1.05, 6e-3)
mv = pd.read_csv(R + 'stage1/mvmr_results.tsv', sep='\t')
direct = mv[(mv.model == 'biot+psa->prac') & (mv.exposure == 'biot')].iloc[0]
check('F4a', 'Direct effect conditional on PSA, OR', np.exp(direct.beta), 1.04, 6e-3)
psa = mv[(mv.model == 'biot+psa->prac') & (mv.exposure == 'psa')].iloc[0]
a = mr('stage1/mr_results.tsv', 'biot->psa', 'IVW').beta
check('F4b', 'PSA channel share, product method', a * psa.beta / mr('stage1/mr_results.tsv', 'biot->prac', 'IVW').beta, 0.50, 6e-3)
co = open(R + 'stage1/coloc_results.txt').read()
import re
pp = [float(x) for x in re.search(r'FKBP4.*?\[([^\]]+)\]', co).group(1).split()]
check('F5', 'FKBP4 colocalization, posterior probability H4', pp[4], 1.0, 0.01)
pp2 = [float(x) for x in re.search(r'2p22.*?\[([^\]]+)\]', co).group(1).split()]
check('F5b', '2p22 colocalization, posterior probability H4', pp2[4], 0.98, 0.01)

# ---------------- Stage 2: rival explanations ----------------
rv = pd.read_csv(R + 'stage2/reverse_results.tsv', sep='\t')
r0 = rv[(rv.analysis == 'prostate cancer liability->biot') & (rv.method == 'IVW')].iloc[0]
check('F6', 'Reverse MR, SD per doubling of odds', r0.beta * np.log(2), -0.0004, 6e-4)
mm = pd.read_csv(R + 'stage2/mvmr_metabolic.tsv', sep='\t')
full = mm[(mm.outcome == 'prac') & (mm.model == 'biot+bmi+igf+fi') & (mm.exposure == 'biot')].iloc[0]
check('F7', 'Metabolic adjusted OR, all three', full.OR, 1.28, 6e-3)
dp = pd.read_csv(R + 'stage2/cml_dp.tsv', sep='\t')
check('F8a', 'MR cML, bioavailable T on PSA', dp.iloc[0].beta, 0.178, 6e-4)
check('F8b', 'MR cML, PRACTICAL, OR', np.exp(dp.iloc[1].beta), 1.16, 6e-3)
up = pd.read_csv(R + 'stage2/uni_metabolic_psa.tsv', sep='\t')
check('F9a', 'BMI on PSA, IVW', up[(up.analysis == 'bmi->psa') & (up.method == 'IVW')].iloc[0].beta, -0.094, 6e-4)
check('F9b', 'IGF-1 on PSA, IVW', up[(up.analysis == 'igf->psa') & (up.method == 'IVW')].iloc[0].beta, 0.003, 6e-4)

# ---------------- Stage 3: synthesis, recomputed from trial counts ----------------
T = pd.read_csv(R + 'stage3/trial_counts_from_sharifan_supp.tsv', sep='\t')
def lc(n, k): return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
g = np.linspace(np.log(0.05), np.log(20), 4001); ll = np.zeros_like(g)
for _, t in T.iterrows():
    m = t.e1 + t.e0; xs = np.arange(max(0, m - t.n0), min(m, t.n1) + 1)
    ll += lc(t.n1, t.e1) + lc(t.n0, t.e0) + t.e1 * g - logsumexp((lc(t.n1, xs) + lc(t.n0, m - xs))[None, :] + np.outer(g, xs), axis=1)
ll -= ll.max(); ci = g[ll >= -1.92]
check('F10', 'Exact conditional OR, 10 trials', np.exp(g[ll.argmax()]), 0.89, 6e-3)
def post(mu0, s0):
    lp = ll + stats.norm.logpdf(g, mu0, s0); p = np.exp(lp - lp.max()); p /= p.sum()
    return np.exp(g[np.searchsorted(np.cumsum(p), 0.5)]), p[g > 0].sum()
sd = (ci.max() - ci.min()) / 3.92
lo, hi = -2.0, 5.0
for _ in range(60):
    mid = (lo + hi) / 2
    lo, hi = (lo, mid) if post(mid, sd)[1] >= 0.95 else (mid, hi)
check('F11', 'Reverse Bayes threshold, prior as firm as trials, OR', np.exp(hi), 2.20, 0.02)
check('F12', 'Posterior OR, PSA channel removed, 1 SD', post(0.0397, 0.1139)[0], 1.02, 6e-3)
P = pd.read_csv(R + 'stage4/pooled_priors.tsv', sep='\t')
tot = P[P.prior == 'total pooled'].iloc[0]; dr = P[P.prior == 'direct pooled'].iloc[0]
check('F14', 'Pooled genetic prior, total, OR', np.exp(tot.beta), 1.149, 6e-3)
check('F15', 'Pooled genetic prior, direct, OR', np.exp(dr.beta), 1.002, 6e-3)
check('F16', 'Posterior, PSA channel removed, pooled prior', post(dr.beta, dr.se)[0], 0.99, 6e-3)
check('F17', 'Posterior, face value, pooled prior', post(tot.beta, tot.se)[0], 1.14, 6e-3)
C3 = pd.read_csv(R + 'stage4/c3_psa_loci_sensitivity.tsv', sep='\t')
e = C3[(C3.outcome == 'prac') & (C3.set == 'excluding PSA biology and PCa loci') & (C3.exposure == 'bioavailable testosterone (direct)')].iloc[0]
check('F18', 'Direct effect excluding PSA biology and PCa loci, OR', e.OR, 1.057, 6e-3)
D5 = pd.read_csv(R + 'stage4/cml_dp_500.tsv', sep='\t')
check('F19', 'MR cML, 500 perturbation replicates, PSA', D5.iloc[0].beta, 0.178, 6e-4)
CA = pd.read_csv(R + 'stage6/cause_results.tsv', sep='\t')
cp = CA[CA.analysis == 'biot->psa'].iloc[0]
check('F20', 'CAUSE, PSA: causal preferred over sharing, delta ELPD', cp.delta_elpd, -5.423, 6e-3)
check('F21', 'CAUSE, PSA: one-sided p', cp.p_one_sided, 0.00087, 6e-5)
check('F22', 'CAUSE, PSA: causal gamma', cp.gamma, 0.15, 6e-3)
cq = CA[CA.analysis == 'biot->prac'].iloc[0]
check('F23', 'CAUSE, prostate cancer: delta ELPD', cq.delta_elpd, -3.982, 6e-3)
check('F24', 'CAUSE, prostate cancer: one-sided p', cq.p_one_sided, 0.0274, 6e-4)
z = (0.1762 - np.log(0.85)) / np.sqrt(0.0511**2 + ((np.log(1.62) - np.log(0.44)) / 3.92)**2)
check('F13', 'Face value genetics versus trial only synthesis, conflict p', 2 * stats.norm.sf(z), 0.314, 2e-3)

w = max(len(r[1]) for r in rows)
for r in rows: print(f'{r[0]:<5} {r[1]:<{w}}  got {r[2]:>9.4f}  expected {r[3]:>8.4f}  {r[4]}')
if fails:
    print('DRIFT in ' + ', '.join(fails), file=sys.stderr); sys.exit(1)
print(f'All {len(rows)} asserted values re-derived within tolerance.')
