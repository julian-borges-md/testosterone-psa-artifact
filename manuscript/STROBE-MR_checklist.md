# STROBE-MR reporting checklist

**Title.** Testosterone raises PSA, and much of its apparent prostate cancer risk runs through that channel: a Mendelian randomization and Bayesian synthesis study

**Author.** Julian Borges, MD, MS

| # | Item | Where reported |
|---|---|---|
| 1 | Title and abstract: state the design | Title names Mendelian randomization and Bayesian synthesis; structured abstract under Background, Methods, Results and Conclusions |
| 2 | Background: why Mendelian randomization is helpful here | Introduction. The androgen saturation model makes a directional prediction that is open to genetic test and that observational data cannot isolate from confounding by adiposity and by screening behaviour. The design bounds rather than locates any threshold, which is stated in the Results and in Table 6 |
| 3 | Objectives: state prespecified causal hypotheses | Methods, Design. One primary test is declared, the effect of bioavailable testosterone on PSA. The saturation decision rule was fixed in advance and encodes sign, significance and estimator agreement only |
| 4 | Study design and data sources | Methods, Design and Data; Table 1 |
| 5 | Setting and participants of contributing studies | Table 1, with sample sizes, ancestry, accessions and collection period. UK Biobank 2006 to 2010, 178,782 to 194,453 men; PSA in 392,522 men; PRACTICAL 79,148 cases; FinnGen release 12, 20,368 cases |
| 6 | Describe the genetic variants used as instruments | Methods, Data; Supplementary Methods section 1, with thresholds, clumping parameters and counts per exposure. Harmonised variant level inputs for all eight exposure to outcome pairs are committed at results/stage7_matrix/harmonised/ |
| 7 | Exposure, outcome and other variables, and how measured | Table 1. Exposures are inverse normalised per standard deviation; PSA is on the natural log scale; cancer outcomes are case control odds ratios |
| 8 | Assess the assumption of association with the exposure | Results, saturation section: mean F 94.1 to 108.2 across the three PSA exposures. Sanderson conditional F of 20.44 and 21.23 for the multivariable models, and 25.9 for SHBG conditional on bioavailable testosterone |
| 9 | Assess independence from confounders | Multivariable MR adjusting for BMI, IGF-1 and fasting insulin; Results, Rival explanations excluded. Reverse causation tested by running genetic liability to prostate cancer against each hormone |
| 10 | Assess the exclusion restriction, no horizontal pleiotropy | Methods, Estimators; Results, Rival explanations excluded. MR Egger intercept 0.0050, p = 1.3 x 10 to the minus 4 for bioavailable testosterone on PSA, indicating directional pleiotropy; CAUSE preferred the causal model over correlated sharing, one sided p = 0.00087, with eta minus 0.30 acting against the causal effect; colocalization at FKBP4 and near SRD5A2. Six of the seven estimators give a positive point estimate, MR Egger alone excepted |
| 11 | Statistical methods, including the estimators used | Methods, Estimators, with the inverse variance weighted estimator stated and explained in narrative; Supplementary Methods section 2 defines all seven. Every estimator is cross validated against a published implementation, Supplementary Table S5 |
| 12 | Assessment of effect modification and interactions | Not performed. Disclosed in Table 6 and in Supplementary Table S9: testing modification by polygenic risk or ancestry requires individual level data, which this study does not use |
| 13 | Sensitivity analyses and additional analyses | Methods, Decomposition and Synthesis; Results, saturation bound, PSA instrument sensitivity under Steiger filtering, the collider negative control, and the randomized trial synthesis. Supplementary Tables S3, S4, S6, S7 and S8 carry the threshold bound on both scales, the SHBG sensitivity of the conversion, the MR RAPS loss sensitivity for every exposure and the contrast, both decomposition methods, and all compatibility specifications live and withdrawn |
| 14 | Software, code and data availability | Methods, Reproducibility; Declarations, Data availability. Three assertion modules run on every push; the full matrix run and the cross validation run on dispatch and schedule. Supplementary Methods section 7 |
| 15 | Descriptive data on the genetic variants | Table 2 and Table 4 give instruments analysed, Cochran Q, I squared and the MR RAPS overdispersion parameter per analysis. Full instrument lists and harmonised inputs are in the repository |
| 16 | Main results with measures of uncertainty | Tables 2, 3, 4 and 5, and Figures 1 and 2; every estimate with a 95 percent interval. Supplementary Tables S1 to S10 carry the reverse Bayes, threshold, sensitivity, cross validation, decomposition, compatibility, limitations and CAUSE results |
| 17 | Results of additional or sensitivity analyses | Results, Rival explanations excluded, Randomized trials show the same channel, and Lethal disease cannot be tested; Supplementary Tables S3, S4, S6, S7, S8 and S10 |
| 18 | Key results in the context of the stated hypothesis | Discussion, opening two paragraphs. The prostate responds to testosterone across the sampled range, and once that PSA response is accounted for the apparent effect on prostate cancer is no longer distinguishable from the null though it is not shown to be zero |
| 19 | Limitations, including assumption violations and the direction of bias | Discussion, Limitations; main Table 6 with 10 rows and the complete table at Supplementary Table S9 with 26 rows, each stating the expected direction of bias |
| 20 | Interpretation, generalisability and funding | Discussion and Conclusion; Declarations, Funding, Conflict of interest and Ethical approval |

## Items requiring specific comment

**Item 2, what the design can and cannot establish.** The manuscript does not claim to falsify the androgen saturation model and states the reason in its own limitations: a linear instrumental variable average over a population in which some men lie below the threshold is compatible with strict saturation. The analysis bounds that possibility, requiring the below threshold stratum to carry a slope about 11.5 times the observed average at 8 nmol/L and about 15,000 times at 4 nmol/L, and a bound is not a rejection. The checklist entry above and the manuscript carry one formulation of this.

**Item 10, exclusion restriction.** The MR Egger intercept is positive and significant for bioavailable testosterone on PSA, which indicates directional pleiotropy and is disclosed rather than corrected after the fact. The manuscript proposes tissue sensitivity loci as a hypothesis and concedes that a reader may regard the negative MR Egger slope as unexplained. Six of seven estimators give a positive point estimate and six of seven intervals exclude the null under the adopted MR RAPS setting, five of seven under an L2 loss.

**Item 11, estimator provenance.** The MR RAPS row is taken from the published implementation, mr.raps 0.2 with a Huber loss and overdispersion enabled, rather than from this study's own code, which was never validated against a reference. MR RAPS is setting dependent on these data and every count that rests on that row states the setting it depends on; Supplementary Tables S6a to S6d.

**Item 12, effect modification.** Not assessed, and not assessable without individual level data. Named as a limitation and as a future direction.

**Item 13, the negative control.** One analysis beyond the standard set is reported: a negative control distinguishing mediation from collider bias when conditioning on PSA. It favours mediation but does not exclude a collider component, so the PSA channel share is reported as a range with an upper bound rather than as a point estimate.

**Item 14, reproducibility.** The repository carries the harmonised variant level inputs, so every cell of the estimator matrices is recomputable offline without access to the source GWAS. Every correction made during development was recorded as it was made, including those identified by the author rather than by a reviewer.

**Item 19, prespecification.** Hypotheses, priors and decision rules were written and time stamped in the committed analysis log before any corresponding estimate was computed. This is internal prespecification and not external registration, and it is disclosed as a limitation in the manuscript rather than described as registration.
