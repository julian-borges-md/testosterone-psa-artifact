# Phase 6 prespecification (CAUSE)
Logged: 2026-10-12T01:54:07Z before running CAUSE.
Purpose: close limitation L3b. CAUSE compares a sharing model (correlated pleiotropy only) against a causal model, and reports the difference in expected log pointwise posterior density (delta ELPD) with a one-sided p value.
Prespecified reading, fixed before computation:
- If the causal model is preferred over the sharing model (p < 0.05) for bioavailable testosterone on PSA, correlated pleiotropy is not a sufficient explanation of the saturation finding, and L3b is closed.
- If the causal model is NOT preferred for bioavailable testosterone on PRACTICAL prostate cancer, that is consistent with the detection reading, in which the cancer association is not a direct causal effect; it is reported as such and not as a failure.
- Any result is reported whatever its direction. Posteriors for the sharing parameter eta and the causal parameter gamma are reported with credible intervals.
Parameters: 1,000,000 variant subset for nuisance parameter estimation (seed 20261012), LD pruning at r2 < 0.01 within 100 kb using the 1000 Genomes EUR panel, p value threshold 1e-3 for pruning, CAUSE 1.2.0 defaults otherwise.

## Result (recorded after computation)
biot->PSA: causal model preferred over sharing, delta ELPD -5.42 (SE 1.73), z = -3.13, p = 0.00087. gamma 0.15 (0.12 to 0.18). Prespecified criterion met: correlated pleiotropy is not a sufficient explanation; L3b closed.
biot->PRACTICAL prostate cancer: causal model preferred, delta ELPD -3.98 (SE 2.07), z = -1.92, p = 0.027. gamma 0.12 (0.05 to 0.19). Weaker evidence than for PSA, and the sharing model is not rejected against the null (z = -1.61).
Software note: CAUSE 1.2.0 required two source patches to run under R 4.3.3 with loo 2.10.1. Both are recorded in docs/CAUSE_PATCHES.md and the patched source is not redistributed.
