**Supplementary material**

**Testosterone Raises PSA, and Much of Its Apparent Prostate Cancer Risk Runs Through That Channel**

Julian Borges, MD, MS.

Every value in this document is re-derived from the committed inputs by the verification modules described in Supplementary Methods, section 7. The repository is at github.com/julian-borges-md/testosterone-psa-artifact, at tag v6.4, archived at https://doi.org/10.5281/zenodo.23352541.

**Supplementary Methods**

**1. Instruments.** Variants at p below 5 times 10 to the minus 8 were restricted to biallelic SNPs with minor allele frequency of at least 1 percent and imputation quality of at least 0.8, then clumped in PLINK 1.90b7.2 against 1000 Genomes European samples (n = 503) at r squared below 0.001 within 10 Mb. This gave 93 bioavailable testosterone, 241 SHBG, 179 total testosterone, 329 PSA, 137 prostate cancer, 78 BMI, 353 IGF-1 and 38 fasting insulin instruments before harmonisation. Harmonisation removed allele mismatches, palindromic SNPs with minor allele frequency above 0.42, and SNPs whose allele frequencies differed by more than 0.2 between sources. Counts at every step are in the committed harmonisation log. Mean F statistics were 94.1 for bioavailable testosterone, 94.6 for total testosterone and 108.2 for SHBG.

**2. Estimator definitions.** The primary estimate is multiplicative random effects IVW. MR Egger provides an intercept test for directional pleiotropy. The weighted median is consistent if at least half the weight comes from valid instruments. The weighted mode uses the zero modal pleiotropy assumption. The Q outlier variant removes variants whose Cochran Q contribution exceeds a Bonferroni threshold and refits IVW. MR cML uses BIC model selection with data perturbation standard errors over 500 replicates. MR RAPS uses a Huber loss with the overdispersion parameter estimated jointly; following the cross validation in section 6 the reported MR RAPS row is taken from the published package mr.raps 0.2 rather than from the study's own implementation.

**3. Steiger filtering.** Applied on the observed scale. The justification differs by analysis and is stated per analysis rather than once: for analyses with PRACTICAL as the outcome the case fraction of 0.564 exceeds population prevalence by more than an order of magnitude, so a liability scale correction would mislead and the observed scale comparison is used; for the PSA analyses the outcome is continuous and the scale question does not arise.

**4. The threshold conversion and its constants.** The saturation threshold of 4 to 8 nmol/L is a serum total testosterone threshold. A lognormal was fitted to a published UK Biobank male total testosterone distribution (mean 12.2, SD 3.5 nmol/L, n = 16,237), giving mu 2.46189 and sigma 0.28123. Conversion to the bioavailable scale used the Vermeulen equation with the association constant for SHBG binding Kt = 1.0 times 10 to the 9 per mol/L, the association constant for albumin binding Ka = 3.6 times 10 to the 4 per mol/L, and an albumin molecular weight of 66,500 g/mol, at the median SHBG of that sample (38.9 nmol/L) and an albumin concentration of 43 g/L. Writing N = 1 + Ka times the albumin molarity, free testosterone solves the quadratic with a = N times Kt and b = N + Kt times (SHBG minus total testosterone), and bioavailable testosterone is N times the free concentration. In words, the equation partitions total testosterone into the fraction bound tightly to SHBG, the fraction bound loosely to albumin and the free fraction, and treats the loosely bound and free fractions together as bioavailable. The calibration check returns a free testosterone of 216.9 pmol/L at the male mean against a published 214.4 pmol/L, a discrepancy of 1.2 percent.

**5. Power.** The noncentrality parameter is lambda = beta times the square root of N times R squared times K times one minus K, where beta is the log odds ratio per SD, N the total sample, R squared the exposure variance explained by the instruments (0.040) and K the case fraction. In words, power rises with sample size, instrument strength and the balance of cases to controls. Scenarios hold the control count between 60,789 and 61,185, the PRACTICAL control count, and vary cases from 1,500 to 79,148, giving case fractions from 0.024 to 0.564. Alpha is 0.05 two sided and equivalence is assessed by two one sided tests against a margin of OR 1.10.

**6. Cross validation against published implementations.** Every estimator was run through a published package on the identical committed harmonised inputs, and the loss sensitivity covers all three PSA exposures and the contrast: MendelianRandomization 0.10.0 for IVW, MR Egger, the weighted median, the weighted mode and MR cML, and mr.raps 0.2 for MR RAPS. Fifty six cells were compared. The closed form estimators agree to machine precision: the maximum absolute difference is 4.4e-16 for IVW, 3.3e-16 for MR Egger, 1.1e-16 for the Q outlier variant and 1.1e-15 for the weighted median, across all eight exposure to outcome pairs. The weighted mode and MR cML differ by at most 0.020 and 0.033 respectively, which is Monte Carlo noise from differing replicate counts, the reference cML using 100 perturbations against this study's 500. MR RAPS differs systematically: the study's own implementation returned a smaller overdispersion parameter than mr.raps in every one of the eight analyses, by 2.58 to 11.11 fold, and on the primary analysis gave 0.1392 with tau squared 1.80e-05 against the published 0.0876 with tau squared 2.00e-04. The direction is systematic rather than random across all eight analyses. The published implementation was adopted for the MR RAPS row throughout on the principle that mr.raps is the estimator this manuscript claims to use and the study's reimplementation was never validated against a reference. The discrepancy is not characterised as a defect in the study code: the published package's own range across loss functions, 0.0089 to 0.1763 on the primary analysis, contains the study value of 0.1392, and the published Tukey setting returns an overdispersion parameter of 1.05e-05, smaller than the study value of 1.80e-05. The defensible statements are that the study implementation was never validated and that MR RAPS is setting dependent on these data (Table S6). The study's own values are retained in the committed matrix and are still recomputed by the assertion module so that the discrepancy remains auditable.

**7. Verification.** Four modules run in continuous integration on every push, without network access. The matrix module performs 166 assertions, of which 141 are recomputations from the committed harmonised variant level data and 25 are transcription checks against committed tables, and it additionally asserts that the overdispersion parameter is strictly positive and that the MR RAPS standard error differs from the fixed effect IVW standard error by more than five percent, which is a regression test for a defect found and corrected during development. The study level module performs 32 assertions, of which 7 are recomputations from trial counts and 25 are transcription checks. A third module parses the manuscript and asserts that all 71 printed matrix values, the bold markings and the estimator counts stated in the text match the committed matrix, that every value plotted in the display figures equals its committed matrix cell and that no figure labels an estimator with a method qualifier the manuscript does not use, that every supplementary table is cited at least once from the main text, that no sentence in the abstract, the Key Messages or the Discussion opening asserts a falsification the limitations table denies, and that every count resting on the MR RAPS row carries the setting it depends on. That last assertion is computed rather than listed: the MR RAPS row is substituted under each loss from the committed sensitivity tables, a statement requires qualification exactly when its count moves, and the alternative must appear within 300 characters of the phrase that prints the statement. The module rejects each earlier draft that carried the defect and passes the present one. A fourth mode of that module, run before submission rather than on every push, additionally asserts that no unresolved placeholder remains in the manuscript. The cross validation, including the loss sensitivity over every exposure and the contrast, runs as a scheduled and dispatchable job because it requires an R toolchain, as does the full matrix run, which additionally recomputes the 300 replicate MR RAPS bootstrap and the 500 replicate MR cML data perturbation standard errors for a total of 182 checks. The bootstrap replay draws its own index stream and is therefore an independent bootstrap at the same seed rather than a bit for bit reproduction, with a Monte Carlo agreement band of 5e-3 and a largest observed residual of 0.0039; the data perturbation replay consumes its generator identically and reproduces to within 1e-6.

**8. Protocol deviations.** Two. CAUSE nuisance parameters were estimated from a random subset of 100,000 variants rather than the prespecified 1,000,000 under a memory constraint, which affects their precision and not the model comparison. The Bayesian prior was changed from the PRACTICAL only estimate to the fixed effect meta-analysis of PRACTICAL and FinnGen, because selecting the stronger of two available estimates is not defensible even when, as here, the stronger estimate works against the authors' own conclusion.

**Supplementary Tables**

**Table S1. Reverse Bayes: minimum prior median odds ratio for the claim to remain at least 95 percent probable.** Firmness is the ratio of the prior standard deviation to the trial likelihood standard deviation of 0.2885 on the log odds scale.

| Prior firmness | Prior SD on the log odds scale | Minimum prior median OR |
|---|---|---|
| Twice as firm as the trials | 0.144 | 1.34 |
| As firm as the trials | 0.288 | 2.20 |
| Half as firm as the trials | 0.577 | 13.3 |

**Table S2. Implications, with the evidence each rests on and the horizon over which it holds.**

| Audience | Implication | Evidence base and horizon |
|---|---|---|
| Clinicians prescribing testosterone | Genetic evidence of risk should not be read as evidence that therapy causes clinically significant cancer; randomized outcomes remain the relevant evidence | Randomized trials in men with symptomatic hypogonadism and a PSA below about 3.0; median follow up about 2 years; does not transport to men near a biopsy threshold |
| PSA monitoring on therapy | A modest rise is expected pharmacodynamics; interpret against a treatment specific baseline rather than population thresholds alone | Genetic bioassay estimate plus trial biopsy data; no direct trial of a monitoring strategy, so this is mechanism based guidance |
| Screening and risk models | PSA based tools and polygenic scores should account for androgen status and BMI, which shift PSA independently of tumour biology | Cross exposure MR within this study; descriptive, three exposures only |
| Biomarker choice | Bioavailable or free testosterone tracks prostate androgen action more closely than total testosterone; the difference is one of magnitude, not of presence against absence | Contrast in Table 3, all seven estimators; European ancestry exposure samples |
| MR methodology | Any MR of an exposure that moves PSA on prostate cancer should report the PSA effect, the screening intensity of the outcome cohort, and aggressive disease results together | Methodological, follows from the decomposition and cohort contrast here |
| Research priority | Aggressive and fatal disease, estimated as a ratio of effects rather than an equivalence test, and functional work at FKBP4 | Requires cohorts not currently public; about 40,000 fatal cases for an equivalence test |

**Table S3. The saturation bound on both scales.** Lognormal fitted to the published UK Biobank male total testosterone distribution. Bioavailable equivalents use the Vermeulen conversion at the median SHBG of 38.9 nmol/L and albumin of 43 g/L, with constants given in Supplementary Methods section 4.

| Total testosterone threshold | Bioavailable equivalent | Fraction of men below | Required slope ratio |
|---|---|---|---|
| 4 nmol/L | 1.60 nmol/L | 0.0065 percent | about 15,000 |
| 6 nmol/L | 2.44 nmol/L | 0.86 percent | about 116 |
| 8 nmol/L | 3.32 nmol/L | 8.69 percent | about 11.5 |

**Table S4. Sensitivity of the bioavailable equivalent of the 8 nmol/L threshold to SHBG.** The mapping between the two scales is not monotone across men, so the equality of tail fractions holds only at fixed SHBG.

| SHBG | Bioavailable equivalent of 8 nmol/L total testosterone |
|---|---|
| 20 nmol/L | 4.74 nmol/L |
| 30 nmol/L | 3.87 nmol/L |
| 38.9 nmol/L (median) | 3.32 nmol/L |
| 50 nmol/L | 2.81 nmol/L |
| 70 nmol/L | 2.20 nmol/L |

**Table S5. Cross validation against published implementations.** Point estimates on identical committed inputs. MendelianRandomization 0.10.0 for all estimators except MR RAPS, which uses mr.raps 0.2.

| Analysis | Estimator | This study | Published implementation | Absolute difference |
|---|---|---|---|---|
| Bioavailable T on PSA | IVW | 0.08478 | 0.08478 | 8.33e-17 |
| Bioavailable T on PSA | MR Egger | -0.03588 | -0.03588 | 2.29e-16 |
| Bioavailable T on PSA | Weighted median | 0.12979 | 0.12979 | 2.50e-16 |
| Bioavailable T on PSA | Weighted mode | 0.14898 | 0.14849 | 4.89e-04 |
| Bioavailable T on PSA | IVW, Q outliers removed | 0.14468 | 0.14468 | 1.11e-16 |
| Bioavailable T on PSA | MR cML BIC | 0.17802 | 0.17382 | 4.20e-03 |
| Bioavailable T on PSA | MR RAPS Huber | 0.13920 | 0.08765 | 5.15e-02 |
| Total T on PSA | IVW | 0.01044 | 0.01044 | 8.15e-17 |
| Total T on PSA | MR Egger | -0.07074 | -0.07074 | 5.55e-17 |
| Total T on PSA | Weighted median | 0.01127 | 0.01127 | 6.42e-17 |
| Total T on PSA | Weighted mode | 0.00963 | 0.01156 | 1.93e-03 |
| Total T on PSA | IVW, Q outliers removed | 0.03803 | 0.03803 | 2.08e-17 |
| Total T on PSA | MR cML BIC | 0.05692 | 0.05747 | 5.55e-04 |
| Total T on PSA | MR RAPS Huber | 0.02157 | -0.01997 | 4.15e-02 |
| SHBG on PSA | IVW | -0.07572 | -0.07572 | 8.33e-17 |
| SHBG on PSA | MR Egger | -0.13190 | -0.13190 | 2.22e-16 |
| SHBG on PSA | Weighted median | -0.06938 | -0.06938 | 6.94e-17 |
| SHBG on PSA | Weighted mode | -0.08062 | -0.07900 | 1.62e-03 |
| SHBG on PSA | IVW, Q outliers removed | -0.06699 | -0.06699 | 2.78e-17 |
| SHBG on PSA | MR cML BIC | -0.06005 | -0.05423 | 5.82e-03 |
| SHBG on PSA | MR RAPS Huber | -0.05844 | -0.07066 | 1.22e-02 |
| Bioavailable T on PRACTICAL | IVW | 0.17616 | 0.17616 | 0.00e+00 |
| Bioavailable T on PRACTICAL | MR Egger | 0.09222 | 0.09222 | 2.78e-17 |
| Bioavailable T on PRACTICAL | Weighted median | 0.11915 | 0.11915 | 4.16e-17 |
| Bioavailable T on PRACTICAL | Weighted mode | 0.10267 | 0.09826 | 4.41e-03 |
| Bioavailable T on PRACTICAL | IVW, Q outliers removed | 0.14204 | 0.14204 | 8.33e-17 |
| Bioavailable T on PRACTICAL | MR cML BIC | 0.14412 | 0.14909 | 4.96e-03 |
| Bioavailable T on PRACTICAL | MR RAPS Huber | 0.15325 | 0.14868 | 4.57e-03 |
| Total T on PRACTICAL | IVW | -0.00190 | -0.00190 | 2.34e-17 |
| Total T on PRACTICAL | MR Egger | 0.01965 | 0.01965 | 1.18e-16 |
| Total T on PRACTICAL | Weighted median | 0.05562 | 0.05562 | 3.54e-16 |
| Total T on PRACTICAL | Weighted mode | 0.06108 | 0.05659 | 4.50e-03 |
| Total T on PRACTICAL | IVW, Q outliers removed | 0.07261 | 0.07261 | 8.33e-17 |
| Total T on PRACTICAL | MR cML BIC | 0.08400 | 0.08043 | 3.57e-03 |
| Total T on PRACTICAL | MR RAPS Huber | 0.05968 | 0.06184 | 2.16e-03 |
| SHBG on PRACTICAL | IVW | -0.17640 | -0.17640 | 0.00e+00 |
| SHBG on PRACTICAL | MR Egger | -0.05705 | -0.05705 | 2.98e-16 |
| SHBG on PRACTICAL | Weighted median | -0.12374 | -0.12374 | 1.05e-15 |
| SHBG on PRACTICAL | Weighted mode | -0.13524 | -0.11480 | 2.04e-02 |
| SHBG on PRACTICAL | IVW, Q outliers removed | -0.08078 | -0.08078 | 8.33e-17 |
| SHBG on PRACTICAL | MR cML BIC | -0.04499 | -0.04964 | 4.65e-03 |
| SHBG on PRACTICAL | MR RAPS Huber | -0.05065 | -0.02085 | 2.98e-02 |
| Bioavailable T on FinnGen | IVW | 0.05236 | 0.05236 | 6.94e-18 |
| Bioavailable T on FinnGen | MR Egger | -0.16846 | -0.16846 | 0.00e+00 |
| Bioavailable T on FinnGen | Weighted median | -0.01630 | -0.01630 | 3.47e-18 |
| Bioavailable T on FinnGen | Weighted mode | -0.04111 | -0.03617 | 4.94e-03 |
| Bioavailable T on FinnGen | IVW, Q outliers removed | -0.02615 | -0.02615 | 9.71e-17 |
| Bioavailable T on FinnGen | MR cML BIC | -0.03973 | -0.00701 | 3.27e-02 |
| Bioavailable T on FinnGen | MR RAPS Huber | -0.00396 | -0.03200 | 2.80e-02 |
| PSA on PRACTICAL | IVW | 1.04209 | 1.04209 | 4.44e-16 |
| PSA on PRACTICAL | MR Egger | 0.89443 | 0.89443 | 3.33e-16 |
| PSA on PRACTICAL | Weighted median | 0.83796 | 0.83796 | 8.88e-16 |
| PSA on PRACTICAL | Weighted mode | 0.77834 | 0.79768 | 1.93e-02 |
| PSA on PRACTICAL | IVW, Q outliers removed | 0.89601 | 0.89601 | 0.00e+00 |
| PSA on PRACTICAL | MR cML BIC | 0.85236 | 0.84755 | 4.81e-03 |
| PSA on PRACTICAL | MR RAPS Huber | 0.86232 | 0.86175 | 5.67e-04 |

**Table S6a. MR RAPS loss function sensitivity, bioavailable testosterone on PSA, with and without overdispersion.** Published implementation, mr.raps 0.2. The package emits a warning that the estimated overdispersion is abnormally large for every row with overdispersion enabled, which is the expected behaviour of a robust profile estimator at an I squared of 90 percent. The Huber row with overdispersion enabled is the setting adopted in Tables 2, 3 and 4.

| Loss function | Overdispersion | Estimate | SE | tau squared |
|---|---|---|---|---|
| L2 | enabled | 0.0089 | 0.0373 | 3.64e-04 |
| L2 | disabled | 0.0939 | 0.0067 | not estimated |
| Huber | enabled | 0.0876 | 0.0293 | 2.00e-04 |
| Huber | disabled | 0.1600 | 0.0072 | not estimated |
| Tukey | enabled | 0.1763 | 0.0111 | 1.05e-05 |
| Tukey | disabled | 0.1800 | 0.0073 | not estimated |

**Table S6b. MR RAPS by exposure and loss function, overdispersion enabled.** Published implementation, mr.raps 0.2, on the committed harmonised inputs. Three statements in the main text besides the primary count rest on the same MR RAPS row, so all three PSA exposures are computed under every loss and each count can be read rather than inferred. The MR RAPS estimate for total testosterone changes sign across the range, negative under L2, negative and not significant under Huber, positive and significant under Tukey, which is the clearest single expression of that instability.

| Exposure | Loss | Estimate | SE | 95 percent interval | Excludes the null | tau squared |
|---|---|---|---|---|---|---|
| Bioavailable T on PSA | L2 | 0.0089 | 0.0373 | -0.0642 to 0.0819 | No | 3.64e-04 |
| Bioavailable T on PSA | Huber (adopted) | 0.0876 | 0.0293 | 0.0302 to 0.1451 | Yes | 2.00e-04 |
| Bioavailable T on PSA | Tukey | 0.1763 | 0.0111 | 0.1546 to 0.1980 | Yes | 1.05e-05 |
| Total T on PSA | L2 | -0.0661 | 0.0325 | -0.1297 to -0.0024 | Yes | 4.69e-04 |
| Total T on PSA | Huber (adopted) | -0.0200 | 0.0243 | -0.0677 to 0.0277 | No | 2.38e-04 |
| Total T on PSA | Tukey | 0.0500 | 0.0120 | 0.0266 to 0.0734 | Yes | 4.13e-05 |
| SHBG on PSA | L2 | -0.1647 | 0.0488 | -0.2603 to -0.0691 | Yes | 2.30e-04 |
| SHBG on PSA | Huber (adopted) | -0.0707 | 0.0330 | -0.1354 to -0.0060 | Yes | 8.64e-05 |
| SHBG on PSA | Tukey | -0.0470 | 0.0219 | -0.0899 to -0.0042 | Yes | 2.80e-05 |

**Table S6c. The bioavailable minus total contrast within each MR RAPS loss.** Difference of the two point estimates with the standard errors combined in quadrature, the same construction as the committed contrast table, the two source studies being independent samples. Only the MR RAPS row of Table 3 depends on the loss; the other six estimators are unchanged.

| Loss | Bioavailable | Total | Difference | SE | z | p |
|---|---|---|---|---|---|---|
| L2 | 0.0089 | -0.0661 | +0.0750 | 0.0494 | +1.52 | 0.1293 |
| Huber (adopted) | 0.0876 | -0.0200 | +0.1076 | 0.0381 | +2.83 | 0.0047 |
| Tukey | 0.1763 | 0.0500 | +0.1263 | 0.0163 | +7.75 | <0.0001 |

**Table S6d. Consequences for every count and range printed in the main text.** Recomputed by substituting the MR RAPS row of Table 2 and Table 3 with each loss setting in turn and leaving the other six estimators unchanged. Three statements move and three do not. The bioavailable testosterone count falls to five of seven under L2. The contrast count falls to five of seven under L2, because the MR RAPS contrast is then 0.0750 with a standard error of 0.0494 at p = 0.129. The total testosterone count rises to three of seven under Tukey, and the number of significantly negative estimates for that exposure rises from one to two under L2; this statement might be thought exempt on the ground that the MR RAPS estimate for total testosterone is negative under the adopted setting, which holds for Huber and L2 but not for Tukey. The SHBG count holds at six of seven under all three losses, so that statement is not setting dependent. The largest p among the significant contrasts is 0.0047 under all three losses: under Huber it is the MR RAPS row at 0.004724 and under L2, where that row is no longer significant, it is the IVW row at 0.004652, which prints identically at four decimal places. The dependency therefore exists in the raw value and does not reach the printed claim.

| Statement as printed in the main text | Huber (adopted) | L2 | Tukey |
|---|---|---|---|
| Bioavailable testosterone raises PSA, six of seven estimators | 6 of 7 | 5 of 7 | 6 of 7 |
| Total testosterone significantly positive on PSA, two of seven | 2 of 7 | 2 of 7 | 3 of 7 |
| Total testosterone significantly negative, MR Egger only | 1 of 7 | 2 of 7 | 1 of 7 |
| SHBG lowers PSA, six of seven estimators | 6 of 7 | 6 of 7 | 6 of 7 |
| Bioavailable exceeds total under all seven estimators | 7 of 7 | 7 of 7 | 7 of 7 |
| Contrast intervals excluding zero, six of seven | 6 of 7 | 5 of 7 | 6 of 7 |
| Largest p among the significant contrasts, 0.0047 | 0.0047 | 0.0047 | 0.0047 |

**Table S7. Decomposition of the testosterone to prostate cancer association, both methods.** The product method is the basis of the reported channel share range. The difference method is reported as evidence of non identification rather than as a candidate estimate, because it attributes the entire unexplained residual to mediation.

| Quantity | Value |
|---|---|
| Total effect, IVW | 0.176162 log odds |
| Direct effect, multivariable MR holding PSA fixed | 0.039655 log odds (OR 1.040, 0.832 to 1.301) |
| PSA on cancer, testosterone held fixed | 1.043614 log odds |
| Indirect effect, product of coefficients | 0.088474 log odds |
| Channel share, product method | 50.2 percent |
| Channel share, product method after Steiger filtering of PSA instruments | 32.5 percent |
| **Reported channel share range (product method)** | **32.5 to 50.2 percent** |
| Channel share, difference method | 77.5 percent, excluded from the range |
| Additivity gap, total minus the sum of direct and indirect | 0.048033 log odds, 27.3 percent of the total effect |
| Conditional F, bioavailable testosterone and PSA on PRACTICAL | 20.44 |
| Conditional F, bioavailable testosterone and PSA on FinnGen | 21.23 |
| Case fraction, PRACTICAL | 0.564 |
| Steiger filtering of PSA instruments | 40 of 312 removed, PSA to cancer coefficient down 35.3 percent |

**Table S8. Compatibility of face value genetic estimates with the randomized evidence.** All specifications, live and withdrawn, as committed in results/stage3/compatibility.tsv.

| Trial anchor | Genetic prior | p | Status |
|---|---|---|---|
| 41 trial exact conditional likelihood | face value, 1 SD | 0.272 | live |
| 41 trial exact conditional likelihood | face value, 2 SD | 0.098 | live |
| Zeng randomized trials only, RR 0.85 | PRACTICAL only, 1 SD | 0.314 | live |
| Zeng randomized trials only, RR 0.85 | PRACTICAL only, 2 SD | 0.139 | live |
| Zeng randomized trials only, RR 0.85 | pooled PRACTICAL and FinnGen, 1 SD | 0.369 | live |
| Zeng all studies including cohorts, RR 0.62 | face value, 1 SD | 0.007 | withdrawn: pools non randomized cohorts |
| Zeng all studies including cohorts, RR 0.62 | face value, 2 SD | 0.001 | withdrawn: pools non randomized cohorts |

The live face value range is therefore 0.098 to 0.369.

**Table S9. Complete limitations table with the expected direction of bias.**

| Issue | Consequence and expected direction |
|---|---|
| Lethal disease untested | Equivalence within OR 1.10 would need about 40,000 fatal cases, at which point the two one sided test reaches 81 percent power |
| A linear IV average cannot falsify a threshold model | A positive average is compatible with saturation if enough men lie below threshold; bounded here but threshold location not tested; non linear MR on individual level data is the decisive test |
| The saturation bound rests on stated assumptions | The threshold is published on the total testosterone scale and converted at fixed median SHBG and albumin; SHBG varies between men so the mapping is not monotone across men. The bound also assumes population share weighting while the estimator uses inverse variance weighting |
| The primary test addresses a population average | Effects in the upper normal range could be overstated if the dose response curve flattens |
| Total testosterone is not null | Significantly positive on PSA in 2 of 7 estimators and on prostate cancer in 2 of 7; the free hormone claim is a magnitude contrast, not presence against absence |
| The fully adjusted metabolic estimate exceeds every other cancer estimate | OR 1.28 at a minimum conditional F of 5.26 admits weak instrument bias, but also constrains any claim of absent effect; not adjudicated |
| Directional pleiotropy and heterogeneity | Egger intercept p = 1.3e-4; I squared 87.9 to 90.0 percent in the PSA analyses; traced in part to tissue sensitivity loci, biasing IVW toward the null and the Egger slope below zero |
| The prespecified rule contained no condition on the Egger intercept or on heterogeneity | Disclosed rather than corrected after the fact |
| The mediated proportion is not identified | Odds ratio non collapsible at a case fraction of 0.564; product and difference methods differ by a factor of 1.54; 27.3 percent of the total effect unaccounted for by their sum. The multivariable direct effect is primary; the channel share is a product method range and an upper bound |
| The collider alternative is not excluded | Two of three negative controls uninformative; no quantitative bias analysis of the tumour to PSA path. The earlier inference that the control ordering favours mediation is withdrawn |
| MR RAPS is unstable on these data, and three printed counts depend on the setting | The published implementation returns 0.009, 0.088 and 0.176 for bioavailable testosterone on PSA under L2, Huber and Tukey losses with overdispersion enabled, and warns that the fitted overdispersion is abnormally large. The bioavailable testosterone count and the contrast count fall to five of seven under L2, and the total testosterone count rises to three of seven under Tukey. The SHBG count and the printed contrast p range hold under all three losses. Tables S6b to S6d |
| The study's own MR RAPS implementation was never validated against a reference | It returned a smaller overdispersion parameter in every analysis, by 2.58 to 11.11 fold; the published implementation is adopted because it is the estimator claimed, and the discrepancy is retained in the committed matrix. The published range across losses contains the study value, so the study value is not characterised as defective |
| PSA GWAS overlaps UK Biobank in its replication stage | Bias toward the observational association, bounded near 0.2 percent of the confounding bias at F = 94 |
| Mixed ancestry across the analysis | European exposures; multiancestry PSA outcome, 28.8 percent non European and 96.5 percent Million Veteran Program; European PRACTICAL; Finnish FinnGen. The exposure to outcome ancestry mismatch attenuates; no transportability claim to African ancestry men |
| Cohort contrast is suggestive | No exposure differed significantly between PRACTICAL and FinnGen, which also differ in power, ancestry and ascertainment |
| Detection signature is descriptive and has an exception | With three exposures the sign concordance test cannot reach significance, and IGF-1 is associated with cancer in FinnGen despite not moving PSA |
| Trial evidence does not transport to men near a biopsy threshold | Contributing trials excluded men with PSA above about 3.0 and men at elevated prostate cancer risk, the stratum where the detection channel would most plausibly operate |
| Event concentration in the trial likelihood | 23 of 50 events from one trial; leave one out moves the likelihood from OR 0.89 to 0.74 but leaves the PSA adjusted posterior at 0.99 |
| Synthesis is conditional on the mediation estimate | Rejecting it returns a posterior of 1.14 and a 12 percent probability of meaningful harm; the saturation and trial biopsy findings are unaffected |
| Exposure is lifelong endogenous testosterone, not administered therapy | Claims about therapy rest on the cited randomized evidence rather than on the genetic estimates |
| Prespecification is internal | Time stamped in the committed analysis log, not registered with OSF or PROSPERO. The PROSPERO number in Table 1 registers the source meta-analysis, not this study |
| Winner's curse | Instruments selected and estimated in the same UK Biobank GWAS; small at mean F above 80 |
| Dose translation | Therapy shift modelled as a 0.5 to 3 SD scenario range rather than estimated; conclusions stable throughout |
| Covariate units | IGF-1 units not stated in the deposited metadata; fasting insulin betas on the raw trait scale |
| CAUSE protocol deviation | Nuisance parameters from 100,000 rather than 1,000,000 variants; affects their precision, not the model comparison |

**Table S10. CAUSE parameters, both models and both outcomes.**

| Analysis | Model | gamma | eta | q |
|---|---|---|---|---|
| Bioavailable T on PSA | sharing | not applicable | 0.19 (0.15 to 0.24) | 0.55 (0.39 to 0.71) |
| Bioavailable T on PSA | causal (preferred, one sided p = 0.00087) | 0.15 (0.12 to 0.18) | -0.30 (-0.41 to -0.19) | 0.08 (0.02 to 0.18) |
| Bioavailable T on PRACTICAL | causal (preferred, one sided p = 0.027) | 0.12 (0.05 to 0.19) | 0.09 | 0.03 |

The negative eta under the preferred causal model for PSA means that whatever correlated pleiotropy remains acts opposite to the causal effect, which would bias the primary estimate toward the null rather than away from it.
**Supplementary Figures**

**Figure S1. Causal model.** Testosterone raises PSA through prostatic androgen receptor signalling, and PSA drives biopsy and diagnosis. The tumour biology path, shown dashed, is not demonstrated for lethal disease.

**Figure S2. Detection signature.** Three exposures against PSA, prostate cancer in PRACTICAL, and prostate cancer in FinnGen. Cancer associations follow each exposure's effect on PSA, with the IGF-1 FinnGen association noted in the main text as an exception.

![](Figure1_conceptual_model.png)

![](Figure3_detection_signature.png)
