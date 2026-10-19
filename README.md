[![verify](https://github.com/julian-borges-md/testosterone-psa-artifact/actions/workflows/verify.yml/badge.svg)](https://github.com/julian-borges-md/testosterone-psa-artifact/actions/workflows/verify.yml)
[![Python](https://img.shields.io/badge/Code-Python-yellowgreen)](#)
[![Design](https://img.shields.io/badge/Design-Mendelian_randomization_%2B_Bayesian_synthesis-blue)](#)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-assertion--checked-brightgreen)](docs/MANIFEST.md)
[![Prespecified](https://img.shields.io/badge/Hypotheses-time_stamped_before_estimation-orange)](prespecification/)
[![License: MIT](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

---

# Testosterone Raises PSA Across the Physiological Range, with Implications for Detection Bias in Prostate Cancer

### Evidence from Mendelian randomization, Bayesian synthesis and randomized trials

**Author:** Julian Borges, MD, MS

---

## Project Overview

For 85 years, clinical practice has carried the belief that testosterone causes prostate cancer. Randomized trials of testosterone therapy do not show it. Mendelian randomization (MR) studies keep finding it. The usual reconciliation, the androgen saturation model, says the prostate stops responding to testosterone above a low threshold, so variation in the normal range cannot matter.

This project tests that reconciliation instead of assuming it. Prostate specific antigen (PSA) is transcribed under direct androgen receptor control, so it serves as an in vivo bioassay of prostatic androgen signaling. If the receptor is saturated, genetically higher testosterone should not move PSA.

**It does.** The prostate keeps responding across the physiological range. That response raises PSA, a raised PSA triggers biopsy, and biopsy finds tumours that were already present. Across exposures, cohorts and trials, the cancer signal follows the PSA signal. Once the PSA channel is removed, the genetic and randomized evidence converge on a posterior centred near the null: the diagnosis arrives earlier rather than the disease being more common. The evidence does not establish absence of effect, and this project does not claim it.

---

## Findings

| # | Finding | Value |
|---|---|---|
| F1 | Bioavailable testosterone raises PSA across the sampled range | **0.085** natural log units per SD (IVW; 95% CI 0.042 to 0.127). Weighted median 0.130, weighted mode 0.149, MR cML 0.178. The prespecified directional criterion is met, and because a linear instrumental variable average over a population in which some men lie below a threshold is compatible with strict saturation, this bounds the saturation model rather than refuting it (see L15) |
| F2 | The prostate reads bioavailable, not total, testosterone | Total testosterone **0.010** (null); SHBG **minus 0.076** |
| F3 | Cancer association present in a screened consortium, not in a registry cohort | PRACTICAL OR **1.19** (1.08 to 1.32); FinnGen OR **1.05** (0.91 to 1.23) |
| F4 | About half of the genetic cancer association runs through PSA | Direct effect conditional on PSA OR **1.04** (0.83 to 1.30); PSA channel share **50%**, to be read as an upper bound (see L12) |
| F5 | Tissue sensitivity loci decouple serum testosterone from prostate signaling | Shared causal variant at **FKBP4** (posterior probability 1.00) and **2p22** near SRD5A2 (0.98) |
| F6 | Reverse causation rejected | Prostate cancer liability on testosterone: **minus 0.0004** SD per doubling of odds |
| F7 | Metabolic confounding rejected | Adjusted for BMI, IGF-1 and fasting insulin: OR **1.28** (no attenuation) |
| F8 | Detection signature across exposures | Exposures that move PSA (testosterone up, BMI down) show cancer effects only in the screened cohort; IGF-1, which does not move PSA, shows the opposite pattern |
| F9 | Trials show the detection channel directly | Therapy **doubles biopsies** (RR 2.38) without increasing cancer (RR 0.85) |
| F10 | The dogma requires an implausible prior | To stay 95% convinced after 41 trials, one must begin believing risk is **doubled** (OR 2.20, exact likelihood) |
| F11 | Synthesis with the PSA channel removed | Posterior OR **0.99** (0.84 to 1.17) using a prior pooled across both cohorts; 1% probability of meaningful harm; stable from 0.5 to 3 SD of therapeutic shift |
| F12b | Correlated pleiotropy excluded by CAUSE | Causal model preferred over the sharing model for PSA (delta ELPD **minus 5.42**, z = minus 3.13, **p = 0.00087**; causal parameter 0.15). For prostate cancer the preference is weaker (p = 0.027), as the detection reading predicts |
| F12 | Robust to PSA instrument choice | Excluding PSA instruments at KLK3, MSMB, TERT, HNF1B and SHBG leaves the direct effect (OR 1.06) and the PSA channel share (0.50) unchanged |

**What this does not show.** It does not show that testosterone has no effect on lethal prostate cancer: no public genetic data set has enough fatal cases to test it (about 40,000 would be needed for an equivalence claim). It does not show that testosterone is safe in men with existing, untreated, castration sensitive cancer, where androgen dependence is established. And the PSA channel is not purely detection, because PSA also indexes prostate androgen activity.

---

## Objectives

1. Test the androgen saturation model causally, using PSA as a bioassay of prostatic androgen receptor signaling.
2. Distinguish bioavailable from total testosterone as the signal the prostate reads.
3. Exclude reverse causation, metabolic confounding and pleiotropy as explanations of the genetic cancer association.
4. Test whether cancer associations across exposures and cohorts are organized by PSA effect, as detection predicts.
5. Quantify how much prior conviction the belief that testosterone therapy causes prostate cancer needs to survive randomized evidence.

---

## Methods & Tools

| Layer | Tool |
|---|---|
| Language | Python 3.11+ (numpy, scipy, pandas, matplotlib); PLINK 1.90b7.2 for LD clumping |
| Instruments | p < 5e-8, biallelic SNPs, MAF of at least 1%, INFO of at least 0.8, r2 < 0.001 within 10 Mb, 1000 Genomes EUR (n = 503) |
| MR estimators | IVW (multiplicative random effects), MR Egger, weighted median, weighted mode, Q outlier removal, MR cML with data perturbation SE, MR RAPS (Huber), Steiger filtering |
| Multivariable MR | Joint clumping, Sanderson conditional F |
| Colocalization | Approximate Bayes factor (single causal variant), 500 kb windows |
| Bayesian synthesis | Exact conditional (noncentral hypergeometric) likelihood across trials; reverse Bayes analysis of credibility |
| Figures | matplotlib, 300 dpi, PNG and vector PDF |

---

## Reproduce

```bash
git clone https://github.com/julian-borges-md/testosterone-psa-artifact.git
cd testosterone-psa-artifact
pip install -r requirements.txt
make verify
```

`make verify` runs the three assertion modules that gate every push. They require no network access and no raw data, and any disagreement with the manuscript exits non-zero.

| Module | What it asserts | Count |
|---|---|---|
| [`src/verify.py`](src/verify.py) | Study level values re-derived from the committed result tables, including the Bayesian synthesis recomputed from the committed trial counts | 32 assertions, 7 recomputed and 25 transcribed |
| [`src/verify_matrix.py`](src/verify_matrix.py) | Every cell of the Tables 2, 3 and 4 estimator matrices, re-derived from the committed harmonised variant level inputs rather than read from a table | 166 checks, 141 recomputed and 25 transcribed |
| [`src/check_manuscript.py`](src/check_manuscript.py) | Every printed matrix value, bold marking, abstract value and estimator count in the manuscript, plus claim consistency, setting qualification, supplementary citation, figure concordance and the journal word limits | 125 manuscript values and gates |

`make selftest` proves the check can fail by perturbing a committed estimate and requiring a non-zero exit. The matrix module additionally carries two regression assertions for the MR RAPS overdispersion defect found in review: that the fitted overdispersion parameter is strictly positive, and that the MR RAPS standard error differs from the fixed effect IVW standard error by more than five percent, which is the signature of a dropped overdispersion term.

Expected output (last lines):

```
All 32 asserted values re-derived within tolerance.
166 checks: 141 recomputed from SNP level inputs, 25 transcribed from committed tables.
All matrix values re-derived within tolerance.
```

`make crossvalidate` runs the comparison against the published packages. It needs an R toolchain, so it runs as a scheduled and dispatchable continuous integration job rather than on every push. Continuous integration status and logs: [Actions](https://github.com/julian-borges-md/testosterone-psa-artifact/actions/workflows/verify.yml).

Full reanalysis from raw summary statistics is described in [`docs/RUNBOOK.md`](docs/RUNBOOK.md). Every claim maps to a command and an expected value in [`docs/MANIFEST.md`](docs/MANIFEST.md).

---

## Repository Contents

```
manuscript/                   the manuscript, supplementary material, title page, cover
                              letter, STROBE-MR checklist and the Word build script
prespecification/             hypotheses and decision rules, time stamped before estimation
src/stage1_saturation/        instruments, harmonization, univariable and multivariable MR, colocalization
src/stage2_detection/         reverse MR, metabolic MVMR, MR cML and MR RAPS
src/stage3_synthesis/         exact likelihood, Bayesian synthesis, reverse Bayes
src/stage4_robustness/        pooled priors, PSA locus sensitivity, cohort heterogeneity, 500-replicate cML
src/stage5_collider/          negative control distinguishing mediation from collider bias
src/stage6_cause/             CAUSE: sharing versus causal model, genome wide
src/stage7_matrix/            regeneration of the full estimator matrix from the harmonised
                              inputs: passA_reproduce.py replays the committed run verbatim,
                              passB_robust.py recomputes the robust estimators with per
                              analysis reseeding, emit_matrices.py writes the display tables
src/stage8_crossvalidation/   comparison against published packages: crossvalidate.R runs
                              MendelianRandomization 0.10.0, raps_reference.R runs mr.raps 0.2,
                              raps_loss_class.R computes the loss sensitivity for every count
                              and range in the manuscript that rests on the MR RAPS row
src/verify.py                 assertion check of every study level manuscript number
src/verify_matrix.py          assertion check of every estimator matrix cell, recomputed from
                              the committed variant level inputs, with the MR RAPS
                              overdispersion regression assertions
src/check_manuscript.py       concordance of the manuscript against the committed matrix, plus
                              the claim consistency, count qualification, supplementary
                              citation and journal limit gates; --submission adds the
                              unresolved placeholder gate
src/selftest.py               proves verify.py fails on perturbed data
results/stage1 .. stage6/     committed result tables for each analysis stage
results/stage7_matrix/        estimator_matrix.tsv with all 56 cells and their diagnostics,
                              contrast_biot_minus_tt.tsv, heterogeneity.tsv, and harmonised/
                              with the variant level inputs for all eight exposure to
                              outcome pairs, so that every cell is recomputable offline
results/stage8_crossvalidation/ crossvalidation.tsv, raps_reference.tsv,
                              raps_loss_sensitivity.tsv, raps_loss_by_exposure.tsv and
                              raps_loss_contrast.tsv
figures/                      display figures (PNG and vector PDF), make_figures.py,
                              make_figure_saturation_full.py, export_figure_data.py and
                              data/ with the exact plotted values
data/download_sources.sh      public sources with checksum verification
docs/                         runbook, manifest, development record, checksums, environment
```

---

## Reproducibility Standard

| Guarantee | Implementation |
|---|---|
| Prespecification | Decision rules written and time stamped before any estimate ([`prespecification/`](prespecification/)); internal, not externally registered |
| Reference verification | Every reference checked against PubMed, the GWAS Catalog, or the publisher of record ([`docs/REFERENCES.md`](docs/REFERENCES.md)); 33 entries, contiguous and in order of first appearance |
| Determinism | Fixed seeds (20261009, 20261011, 7); deterministic exact likelihood for the synthesis |
| Provenance | SHA256 checksums on every input ([`docs/DATA_CHECKSUMS.txt`](docs/DATA_CHECKSUMS.txt)); md5 matched against EBI deposits |
| Claim traceability | Every claim maps to a command and an expected value ([`docs/MANIFEST.md`](docs/MANIFEST.md)) |
| Drift detection | Expected values compiled into `src/verify.py` as assertions; failure exits non-zero |
| Figure auditability | Every plotted value exported to [`figures/data/`](figures/data/); figures redraw from `results/` with no retyped numbers. [`src/check_manuscript.py`](src/check_manuscript.py) asserts that each plotted value equals its committed matrix cell and that no figure labels an estimator with a method qualifier the manuscript does not use after the figures were found to be the one artifact class no module examined |
| Cell level recomputation | Every estimator matrix cell is recomputed from the committed harmonised variant level inputs in [`results/stage7_matrix/harmonised/`](results/stage7_matrix/harmonised/), not transcribed from a result table |
| Independent implementation | Every estimator is cross validated against a published package on the identical committed inputs ([`src/stage8_crossvalidation/`](src/stage8_crossvalidation/)); the MR RAPS row is taken from mr.raps 0.2 rather than from this study's code |
| Regression assertions | The two signatures of the R2 overdispersion defect are asserted on every push: strictly positive fitted overdispersion, and an MR RAPS standard error distinguishable from the fixed effect IVW standard error |
| Claim consistency | [`src/check_manuscript.py`](src/check_manuscript.py) asserts mechanically that no sentence in the abstract, the Key Messages or the Discussion opening asserts a falsification the limitations table denies |
| Setting qualification | The same module asserts that every count resting on the MR RAPS row states the setting it depends on. Which statements require it is computed from the committed loss sensitivity tables rather than listed, so a regeneration that makes a stable count setting dependent starts failing the gate; and the qualification must appear within 300 characters of the phrase that prints the statement, because a region wide search passes by accident. The gate rejects a draft that carried the defect and passes the present one |
| Supplementary citation | Every supplementary table is asserted to be cited at least once from the main text |
| Falsifiability of the check | `src/selftest.py` requires verify.py to fail on perturbed input |
| Error disclosure | development record (``) |

### Known limitations, disclosed

This repository documents its own defects rather than omitting them.

| ID | Limitation |
|---|---|
| L1 | Pipeline scripts use the absolute paths of the original execution environment (`/tmp/gw`, `/tmp/tools`). [`docs/RUNBOOK.md`](docs/RUNBOOK.md) reproduces them with symbolic links. Refactoring to configurable paths is pending |
| L2 | Lethal prostate cancer is untested. Equivalence at OR 1.10 needs about 40,000 fatal cases; no public GWAS approaches this |
| L3 | Prespecification is internal: decision rules are time stamped in the committed log and repository, but not registered with OSF or PROSPERO |
| L3b | CAUSE nuisance parameters were estimated from 100,000 rather than 1,000,000 variants (3 GB memory limit); this affects their precision, not the model comparison. CAUSE 1.2.0 required two source patches to run under R 4.3.3 ([docs/CAUSE_PATCHES.md](docs/CAUSE_PATCHES.md)) |
| L4 | Significant MR Egger intercept for testosterone on PSA (p = 1e-4), traced to tissue sensitivity loci that bias IVW toward the null |
| L5 | The PSA GWAS replication stage includes UK Biobank, overlapping the exposure samples; bias bounded near 0.2% of the confounding bias at F = 94 |
| L6 | Trial level counts were read from the forest plot in a published commentary supplement, available only for the 10 trials with events; they reproduce the published pooled estimate |
| L7 | Translation from lifelong genetic exposure to two years of therapy is modeled as a 0.5 to 3 SD scenario range, not estimated |
| L8 | IGF-1 units are not stated in the deposited metadata; fasting insulin effects use raw trait betas with p values from an inverse normalized analysis |
| L9 | European ancestry exposure GWAS only; no transportability claim to African ancestry men |
| L10 | Winner's curse: instruments were selected and estimated in the same UK Biobank GWAS, inflating apparent instrument strength; small at mean F above 80 |
| L11 | The detection signature is descriptive: with three exposures the sign concordance test cannot reach significance (minimum attainable one-sided p = 0.125) |
| L12 | Conditioning on PSA may induce collider bias, since PSA lies downstream of undiagnosed tumour. A negative control (results/stage5) favours mediation but does not exclude a collider component, so the PSA channel share is an upper bound |
| L13 | The synthesis is conditional on the mediation estimate: rejecting it returns a posterior OR of 1.14 and a 12% probability of meaningful harm. The PSA response to testosterone and the trial biopsy finding do not depend on it |
| L14 | The exposure is lifelong endogenous testosterone, not administered therapy; claims about therapy rest on cited randomized evidence |
| L15 | A linear instrumental variable average cannot refute a threshold model. A positive average is compatible with saturation if enough men lie below the threshold; that possibility is bounded here, requiring the below threshold stratum to carry a slope about 11.5 times the observed average at 8 nmol/L and about 15,000 times at 4 nmol/L, and a bound is not a rejection. Threshold location is not tested; non linear MR on individual level data is the decisive test |

**Corrections are on record.** An apparent conflict between genetic and trial evidence was driven by a pooled estimate that included a registry cohort, and was withdrawn. The genetic prior was later changed from PRACTICAL only to a meta-analysis of both cohorts, which was the more conservative choice against our own hypothesis.

---

## Data Availability

| Source | Access | Identifier |
|---|---|---|
| Bioavailable testosterone, SHBG, total testosterone, men (UK Biobank, Ruth 2020) | Open, GWAS Catalog | GCST90012103, GCST90012109, GCST90012113 |
| PSA levels (Hoffmann 2025) | Open, GWAS Catalog | GCST90461907 |
| Prostate cancer, PRACTICAL OncoArray (Schumacher 2018) | Open, GWAS Catalog | GCST006085 |
| Prostate cancer, FinnGen R12 | Open | C3_PROSTATE_EXALLC |
| BMI (Locke 2015), IGF-1 (Sinnott Armstrong 2021), fasting insulin (Chen 2021) | Open, GWAS Catalog | GCST002783, GCST90019511, GCST90002238 |
| LD reference | Open, MRC IEU | 1000 Genomes phase 3 EUR |
| Randomized trial counts | Published | García Becerra 2026 meta-analysis; Sharifan 2026 supplement; Zeng 2022 |

No controlled access or patient identifiable data is used or redistributed here. Raw summary statistics are not redistributed; `data/download_sources.sh` retrieves and verifies them.

---

## Citation

> Borges, Julian. *Testosterone Raises PSA Across the Physiological Range, with Implications for Detection Bias in Prostate Cancer.* 2026. https://github.com/julian-borges-md/testosterone-psa-artifact

Archived at [https://doi.org/10.5281/zenodo.23483135](https://doi.org/10.5281/zenodo.23483135) (version v6.6). The concept DOI [10.5281/zenodo.23324218](https://doi.org/10.5281/zenodo.23324218) always resolves to the latest version. Machine readable metadata in [`CITATION.cff`](CITATION.cff).

---

## Author

**Julian Borges, MD, MS**
Department of Computer Science, Boston University
ORCID [0009-0001-9929-3135](https://orcid.org/0009-0001-9929-3135)
jyborges@bu.edu

---

## Keywords

`testosterone` • `prostate cancer` • `PSA` • `Mendelian randomization` • `androgen saturation model` • `detection bias` • `Bayesian synthesis` • `reverse Bayes` • `reproducibility` • `endocrinology`

---

## License

Code is licensed under the [MIT License](LICENSE). Figures and documentation are released under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). Upstream datasets retain their original licenses.

---

<div align="center">

**Frontier Translational Research Lab**

Department of Computer Science · Boston University · Harvard Medical School GCSRT Alumni

[![Lab Website](https://img.shields.io/badge/Lab-frontier--lab-002244?style=flat-square)](https://julian-borges-md.github.io/frontier-lab/)
[![BU CS](https://img.shields.io/badge/BU-Computer_Science-cc0000?style=flat-square)](https://www.bu.edu/cs/)
[![HMS](https://img.shields.io/badge/HMS-GCSRT_Alumni-a51c30?style=flat-square)](https://ghsm.hms.harvard.edu/education/global-clinical-scholars-research-training)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0001--9929--3135-a6ce39?style=flat-square&logo=orcid&logoColor=white)](https://orcid.org/0009-0001-9929-3135)
[![CV](https://img.shields.io/badge/Academic_CV-research--profile-4f46e5?style=flat-square)](https://julian-borges-md.github.io/research-profile/)

*Julian Borges, MD, MS · jyborges@bu.edu*

</div>
