# RLO_TPC | Reproducibility Manifest

Every claim maps to a command and an expected value. `make verify` checks all of them in one run.

## Claim to evidence map

| # | Claim | Command | Expected |
|---|---|---|---|
| F1 | Bioavailable testosterone raises PSA across the sampled range, which bounds the saturation model rather than refuting it | `make verify` (F1a to F1c) | IVW 0.085, weighted median 0.130, weighted mode 0.149 |
| F2 | Prostate reads bioavailable, not total, testosterone | `make verify` (F2a, F2b) | Total testosterone 0.010; SHBG minus 0.076 |
| F3 | Cancer association in PRACTICAL, not FinnGen | `make verify` (F3a, F3b) | OR 1.19 and 1.05 |
| F4 | About half of the association runs through PSA | `make verify` (F4a, F4b) | Direct OR 1.04; PSA share 0.50 |
| F5 | Shared causal variants at FKBP4 and 2p22 | `make verify` (F5, F5b) | Posterior probability H4 1.00 and 0.98 |
| F6 | Reverse causation rejected | `make verify` (F6) | Minus 0.0004 SD per doubling of odds |
| F7 | Metabolic confounding rejected | `make verify` (F7) | Fully adjusted OR 1.28 |
| F8 | Robust estimators agree | `make verify` (F8a, F8b) | MR cML 0.178 (PSA), OR 1.16 (PRACTICAL) |
| F9 | Detection signature: BMI lowers PSA, IGF-1 does not move it | `make verify` (F9a, F9b) | Minus 0.094 and 0.003 |
| F10 | Exact pooled trial OR | `make verify` (F10) | 0.89 |
| F11 | Reverse Bayes threshold | `make verify` (F11) | OR 2.20 |
| F12 | Posterior with the PSA channel removed | `make verify` (F12) | 1.02 |
| F13 | No statistical conflict with trial only evidence | `make verify` (F13) | p = 0.314 |
| F14 | Pooled genetic prior across cohorts | `make verify` (F14, F15) | Total OR 1.149; direct OR 1.002 |
| F16 | Posterior with pooled prior, PSA channel removed | `make verify` (F16, F17) | 0.99 and 1.14 |
| F18 | Robust to PSA instrument choice | `make verify` (F18) | Direct OR 1.057 excluding PSA biology and PCa loci |
| F19 | MR cML stable at 500 replicates | `make verify` (F19) | 0.178 |
| F20 | CAUSE: correlated pleiotropy excluded for PSA | `make verify` (F20 to F22) | delta ELPD minus 5.423, p = 0.00087, gamma 0.15 |
| F23 | CAUSE: prostate cancer, weaker preference | `make verify` (F23, F24) | delta ELPD minus 3.982, p = 0.0274 |
| T1 | Trials double biopsies without more cancer | Published input (Zeng 2022, randomized trials) | RR 2.38 and 0.85 |

## Figures

| Figure | Code | Plotted values |
|---|---|---|
| 1 conceptual model | `figures/make_figures.py` | `figures/data/Figure1_conceptual_model_edges.tsv` |
| 2 saturation test | same | `figures/data/Figure2_saturation_test_data.tsv` |
| 3 detection signature | same | `figures/data/Figure3_detection_signature_data.tsv` |
| 4 trials and synthesis | same | `figures/data/Figure4A_trial_inputs.tsv`, `figures/data/Figure4B_posterior_curves.tsv` |

`python3 figures/export_figure_data.py` rewrites the data files from `results/`.

## Verification status

| Category | Status |
|---|---|
| Manuscript values | 32 of 32 re-derived within tolerance by assertion |
| Bayesian synthesis | Recomputed from committed trial counts, not read from stored output |
| Falsifiability of the check | `make selftest` confirms a perturbed estimate triggers DRIFT and a non-zero exit |
| Full pipeline from raw data | Documented in RUNBOOK.md; not run in CI (large downloads, PLINK, LD reference) |
