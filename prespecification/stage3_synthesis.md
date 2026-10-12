# Paper 3 prespecification (RLO_TPC_DogmaCredibility_v1.0)
Logged: 2026-10-11T17:15:01Z, before computation.
Question: given randomized evidence, how strong must a prior belief be for the claim "testosterone therapy increases prostate cancer" to remain credible?
Data (verified via PubMed): Garcia-Becerra 2026 meta-analysis of 41 RCTs (n=11,161): any prostate cancer event OR 0.88 (0.52-1.51); clinically significant OR 1.13 (0.39-3.26). TRAVERSE high grade: 5/2,596 vs 3/2,602.
Method: normal approximation on log OR; conjugate updating; reverse Bayes: minimum prior median OR (at stated prior certainty) needed for posterior P(OR>1) >= 0.95.
Priors: vague N(0,10^2); dogma priors with median OR 1.5 and 2.0; genetics informed priors from Paper 1 MR (OR 1.19 per SD, SE 0.0511 log scale) scaled to 1 and 2 SD shifts.
Thresholds: clinically meaningful harm OR 1.2.
