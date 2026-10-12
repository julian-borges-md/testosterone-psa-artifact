**Testosterone raises PSA, and much of its apparent prostate cancer risk runs through that channel: a Mendelian randomization and Bayesian synthesis study**

**Julian Borges, MD, MS**

Department of Computer Science, Boston University Metropolitan College, Boston, MA 02215, United States

Correspondence: jyborges@bu.edu | ORCID 0009-0001-9929-3135

**Abstract**

**Background.** Genetically higher testosterone predicts more prostate cancer diagnoses, yet randomized trials show no excess. The androgen saturation model reconciles this by holding that prostatic androgen receptors saturate above 4 to 8 nmol/L, so normal range variation cannot affect it; it has never been tested causally. Prostate specific antigen (PSA) is transcribed under androgen receptor control, making it both a bioassay of prostatic androgen signalling and the biopsy trigger.

**Methods.** Two sample Mendelian randomization using male specific instruments for bioavailable testosterone, total testosterone and SHBG from UK Biobank (2006 to 2010; 178,782 to 194,453 men), with PSA in 392,522 men and prostate cancer from PRACTICAL (79,148 cases) and FinnGen (20,368 cases). Seven estimators were computed for every exposure and outcome, followed by multivariable MR and Bayesian synthesis with trial evidence.

**Results.** Bioavailable testosterone raised PSA (inverse variance weighted 0.085 natural log units per SD, 95% CI 0.042 to 0.127), with intervals excluding the null in six of seven estimators (five under an alternative MR RAPS loss). Testosterone raised diagnosed cancer in PRACTICAL (OR 1.19, 1.08 to 1.32) but not FinnGen (OR 1.05, 0.91 to 1.23); conditioning on PSA reduced the direct effect to OR 1.04 (0.83 to 1.30), with 32.5% to 50.2% carried by the PSA channel. In trials, therapy doubled biopsies (RR 2.38, 1.01 to 5.60) without increasing cancer (RR 0.85, 0.44 to 1.62).

**Conclusions.** Testosterone raises PSA, and much of its apparent prostate cancer risk runs through that channel rather than tumour biology. The evidence does not establish absence of effect.

**Keywords:** testosterone; prostate cancer; prostate specific antigen; Mendelian randomization; androgen saturation model; detection bias; Bayesian synthesis.

**Key Messages**

- Genetically higher bioavailable testosterone raises prostate specific antigen across the sampled physiological range, which is to our knowledge the first causal test of the androgen saturation model in humans, although a linear instrumental variable average can bound but not locate the threshold the model posits.
- Between a third and a half of the genetic association between testosterone and diagnosed prostate cancer is carried by the PSA channel, so part of that association reflects when tumours are found rather than whether they arise.
- The residual direct effect is not distinguishable from the null but is not shown to be zero, and an absence claim would require an equivalence test that no available data set can support.

**Introduction**

In 1941 Huggins, Stevens and Hodges showed that castration made metastatic prostate cancer regress [1]. The complementary claim, that giving testosterone accelerates cancer, rested on far weaker evidence and was traced by a later historical review to a single patient measured with an unreliable marker [2]. The belief nonetheless shaped prescribing for decades.

Three bodies of evidence now disagree. Pooled prospective cohorts find no association between endogenous testosterone and prostate cancer across the normal range [3,4]. Randomized trials report no excess cancer [5,6]. Yet Mendelian randomization (MR) studies consistently find that genetically higher free or bioavailable testosterone raises risk, with 29 such studies identified in a 2025 review [7]. The androgen saturation model is the standard reconciliation [8]: prostatic androgen receptors are held to be saturated at roughly 4 to 8 nmol/L of serum total testosterone, so replacement into the normal range would be harmless by mechanism. The model is widely cited in clinical guidance and has never been tested causally.

It makes a directional prediction open to genetic test. PSA is transcribed from KLK3 under direct androgen receptor control, so PSA is an in vivo readout of prostatic androgen signalling: if receptors are saturated, genetic variation in bioavailable testosterone within the physiological range should not move PSA. PSA also links hormone biology to diagnosis: in screened populations most prostate cancers are found when PSA crosses a biopsy threshold. Anything that raises PSA can therefore raise diagnosed cancer without changing tumour biology; obesity illustrates the reverse direction, lowering PSA through plasma hemodilution [9].

We tested the saturation prediction, determined whether the prostate responds to bioavailable or to total testosterone, excluded reverse causation and metabolic confounding, and quantified how much prior conviction the belief that therapy causes prostate cancer requires in order to survive the randomized evidence.

**Methods**

Full methods, including instrument construction, all estimator definitions, the threshold conversion and its constants, the power calculation and the protocol deviations, are in Supplementary Methods.

**Design.** Two sample MR on public summary statistics, followed by Bayesian synthesis with randomized trial evidence. Hypotheses, priors and decision rules were written and time stamped in the committed analysis log before any corresponding estimate was computed. This is internal prespecification, not external registration, and is disclosed as a limitation. One primary test is declared: the effect of bioavailable testosterone on PSA. The saturation decision rule was fixed in advance and encodes sign, significance and estimator agreement only; it contains no condition on the MR Egger intercept or on heterogeneity, and it does not test the location of any threshold.

**Data.** Exposures are the male specific UK Biobank hormone GWAS [10]; the outcome bioassay is a multiancestry PSA GWAS in 392,522 men, of whom 28.8% are non European and 96.5% are drawn from the Million Veteran Program [11]; cancer outcomes are PRACTICAL [12] and FinnGen release 12 [13]; covariates are GIANT BMI [14], UK Biobank IGF-1 [15] and MAGIC fasting insulin [16]; trial evidence comes from a 41 trial meta-analysis [17] with counts from its correspondence supplement [18] and from the randomized trial only stratum of a network meta-analysis [19]. Table 1 gives every source with its ancestry, units and denominators. Instrument counts are 93 bioavailable testosterone, 241 SHBG and 179 total testosterone before harmonisation, and 87, 232 and 173 after it; analysed counts appear in every table.

**Estimators.** The primary estimate is multiplicative random effects inverse variance weighted (IVW) MR [20], which pools the variant specific Wald ratios weighted by precision: each variant contributes how far it moves the outcome divided by how far it moves the exposure, and variants measured more precisely count more in the average. Six further estimators were computed for every exposure and both outcomes, no cell omitted: MR Egger [21], the weighted median [22], the weighted mode [23], IVW after Cochran Q outlier removal, MR cML with data perturbation standard errors from 500 replicates [24], and MR RAPS with a Huber loss and overdispersion [25]. Every estimator was cross validated against a published implementation on the identical committed inputs (Supplementary Table S5): the closed form estimators agree with MendelianRandomization 0.10.0 to machine precision, while our own MR RAPS code, which was never validated against a reference implementation, returned a smaller overdispersion parameter in every analysis, by 2.6 to 11.1 fold, so the reported MR RAPS row is taken from mr.raps 0.2 rather than from our implementation. MR RAPS is setting dependent on these data and the loss function sensitivity is reported in Supplementary Table S6. The fitted tau squared is reported in Tables 2 and 4. CAUSE [26] was run genome wide, and Steiger filtering [27] on the observed scale, justified per analysis.

**Decomposition.** Multivariable MR regressed variant outcome effects jointly on bioavailable testosterone and PSA effects; the coefficient on testosterone is its direct effect holding PSA fixed, and instrument strength used the Sanderson conditional F statistic [28]. The PSA channel was estimated by the product of coefficients. The difference method is reported alongside it, but the two are not equal here and neither identifies a mediated proportion without an additivity assumption on the log odds scale, which the additivity gap shows to be violated at a case fraction of 0.564. The multivariable direct effect, which does not require that assumption, is therefore the primary quantity. Both decompositions are given in Supplementary Table S7.

**Synthesis.** Counts from the 10 trials reporting prostate cancer events (8,631 participants, 50 events) were analysed with an exact conditional noncentral hypergeometric likelihood. The genetic estimate entered as an informative prior, prespecified as the fixed effect meta-analysis of the PRACTICAL and FinnGen estimates. A reverse Bayes analysis [29] then asked how firm a prior belief must be for the claim that therapy increases cancer to remain at least 95% probable after the trials (Supplementary Tables S1 and S2).

**Reproducibility.** Analysis code, committed result tables, the harmonised variant level inputs for all eight exposure to outcome pairs and the verification modules are openly released under an MIT licence at github.com/julian-borges-md/testosterone-psa-artifact, at tag v6.1, archived at https://doi.org/10.5281/zenodo.23324744. Three assertion modules run in continuous integration on every push without network access: 166 matrix assertions, 141 of them recomputed from the committed variant level data, 32 study level values, and a concordance module asserting that every printed, abstract and plotted value, and every supplementary cross reference, matches the committed matrix. The cross validation runs as a scheduled job. Reporting follows STROBE-MR [30].

**Results**

**The prostate responds to testosterone across the sampled range.** Genetically higher bioavailable testosterone raised PSA, with intervals excluding the null in six of seven estimators (Table 2, Figure 1). The sixth is MR RAPS under the Huber loss adopted throughout; under an L2 loss its interval includes the null and the count is five of seven (Supplementary Table S6). The primary IVW estimate of 0.085 corresponds to 8.8% higher PSA per SD. The robust estimators span 0.088 to 0.178, with MR RAPS at the primary IVW estimate rather than above it and MR cML highest. MR Egger is the only estimator not positive. Instruments were strong (mean F 94.1 to 108.2) and heterogeneity high throughout (I squared 87.9 to 90.0%), which is why the estimators disperse and why it is reported per analysis.

A positive average is compatible with strict saturation if enough men lie below the threshold, so we bounded that possibility on the scale the threshold is stated on. Fitting a lognormal to a published UK Biobank male total testosterone distribution (mean 12.2, SD 3.5 nmol/L) [31] gives 0.0065% of men below 4 nmol/L, 0.86% below 6 and 8.69% below 8. For the observed average to arise solely from men below threshold, their slope would need to be about 11.5 times the average at 8 nmol/L and 15,000 times at 4. Applying the Vermeulen equation [32] at the median SHBG of that sample, those thresholds correspond to 1.60, 2.44 and 3.32 nmol/L of bioavailable testosterone; the conversion, its constants and its SHBG sensitivity are in Supplementary Methods and Tables S3 and S4, with a calibration check returning 216.9 pmol/L of free testosterone against a published 214.4. The bound assumes population share weighting of stratum specific slopes whereas the estimator uses inverse variance weighting; the direction survives any reasonable weighting but the multipliers are approximations. This bounds rather than tests the threshold, and distinguishing strict saturation from diminishing returns requires non linear MR on individual level data.

**The prostate responds more strongly to bioavailable than to total testosterone.** Total testosterone is not null: significantly positive on PSA under two of seven estimators, three under a Tukey MR RAPS loss, and significantly negative under MR Egger, and under two of seven estimators at an L2 MR RAPS loss. The defensible statement is a magnitude contrast, tested directly (Table 3). Bioavailable testosterone exceeds total testosterone under all seven estimators, with intervals excluding zero under six, five at an L2 MR RAPS loss, and p from below 0.0001 to 0.0047 under every loss (Supplementary Table S6); under MR Egger both estimates are negative, so the comparison is not between positive effects and the difference is not significant. We do not report the ratio of the two estimates, 3.1 to 15.5 fold, because the IVW denominator of 0.010 (SE 0.015) is indistinguishable from zero. SHBG lowers PSA in six of seven estimators under every MR RAPS loss and retains its effect conditional on bioavailable testosterone (−0.066, −0.123 to −0.010; conditional F 25.9).

The largest heterogeneity outliers paired lower serum testosterone with higher PSA and colocalized with PSA at FKBP4 (posterior probability 1.00), encoding an androgen receptor cochaperone that enhances receptor transactivation without altering hormone binding or nuclear translocation [33], and near SRD5A2 at 2p22 (0.98). These selections were post hoc. CAUSE preferred the causal model over correlated sharing for PSA (one sided p = 0.00087); under that model the correlated pleiotropy parameter eta is −0.30 (−0.41 to −0.19), acting opposite to the causal effect and therefore biasing the primary estimate toward the null (Supplementary Table S10).

**Prostate cancer: a signal that runs through PSA.** Bioavailable testosterone raised prostate cancer risk in PRACTICAL in five of seven estimators and in none in FinnGen (Table 4). Total testosterone is significantly positive on PRACTICAL in two of seven. Conditioning on PSA reduced the direct effect to OR 1.04 (0.83 to 1.30) at a conditional F of 20.4, an interval still compatible with a clinically meaningful increase.

Three values exist for the PSA channel share and are not interchangeable. The product of coefficients gives 50.2%, and the same product after Steiger filtering of the PSA instrument set, which removes 40 of 312 instruments and lowers the PSA to cancer coefficient by 35.3%, gives 32.5%. We report **32.5% to 50.2% as the channel share**, and the term upper bound applies to that range. This range is the product of coefficients across instrument specifications. The difference method gives 77.5%, which lies above it, and is deliberately excluded rather than overlooked: the indirect and direct components sum to 0.128 against a total effect of 0.176, an additivity gap of 0.048 log odds or 27.3% of the total, and that gap discredits the difference method, which attributes the entire unexplained residual to mediation. It is reported as evidence of non identification, not as a candidate value.

The negative control intended to separate mediation from collider bias does not settle the question. Bioavailable testosterone attenuated by 77% on PSA adjustment and IGF-1, which does not move PSA, by 27%, but the IGF-1 univariable estimate is itself null in that instrument set (p = 0.30) and the BMI estimate is too imprecise to interpret, so two of three controls are uninformative. An earlier version read the ordering as favouring mediation; that inference is withdrawn, and a collider component cannot be excluded.

**Rival explanations excluded.** Genetic liability to prostate cancer did not change any hormone. Adjustment for BMI, IGF-1 and fasting insulin did not attenuate the testosterone effect but strengthened it, to OR 1.28 (1.08 to 1.52), the largest and most precisely positive cancer estimate in this study. We report it as such. Two readings are available and we do not adjudicate: the model has a minimum conditional F of 5.26, so weak instrument bias is possible, or the direct effect is not null. Either way it constrains how strongly an absence of effect can be asserted.

**Randomized trials show the same channel.** Therapy more than doubled prostate biopsies (RR 2.38, 1.01 to 5.60) without increasing cancer (RR 0.85, 0.44 to 1.62) (Figure 2A); both come from the randomized trial only stratum of a network meta-analysis whose full network also pools non randomized cohorts [19], and both intervals are wide, so direction is firmer than magnitude. The exact conditional likelihood across the 10 trials with events gave OR 0.89 (0.50 to 1.56). To remain 95% convinced of harm afterwards, a reader must begin believing risk is at least doubled, held as firmly as the trial evidence itself. With the PSA channel removed the posterior is OR 0.99 (0.84 to 1.17), with a 1% probability of meaningful harm (Table 5, Figure 2B); excluding the largest trial, which supplies 23 of 50 events, leaves it at 0.99. Face value genetic estimates are not statistically incompatible with the trials across five specifications (p = 0.098 to 0.369; Supplementary Table S8), so the case against reading them as tumour biology rests on mechanism rather than statistical conflict.

**Lethal disease cannot be tested.** Detecting an OR of 1.19 for aggressive disease needs about 10,000 cases, but equivalence within an OR of 1.10 would need about 40,000, at which point the two one sided test reaches 81% power. No public data set approaches this, so the prespecified equivalence test was not executed.

**Discussion**

Two findings sit at the centre of this study. The prostate responds to testosterone across the sampled range: genetic increases in bioavailable testosterone raise PSA, a direct androgen receptor product, with agreement across six of seven estimators, five under an L2 MR RAPS loss. To our knowledge this is the first causal test of the saturation model in humans, and it bounds rather than locates the threshold. Second, once that PSA response is accounted for, the apparent effect of testosterone on prostate cancer is no longer distinguishable from the null, though it is not shown to be zero.

The dogma holds that testosterone causes prostate cancer. This study does not refute it outright; it shows that much of the genetic evidence for it measures detection rather than disease. Testosterone raises PSA, a raised PSA brings forward biopsy, and biopsy finds a tumour already present, so what changes is the timing of discovery rather than the biology. Three lines support this reading, with differing strength: the mediation estimate, which is an upper bound resting on an assumption shown here to be violated; the pattern across exposures, which is descriptive and has one exception in IGF-1, associated with cancer in FinnGen despite not moving PSA; and the randomized trials, in which therapy doubled biopsies without increasing cancer, which is the most direct of the three.

Morgentaler and Traish reached a similar clinical conclusion by arguing that the prostate stops responding to testosterone [8]. Our evidence says the opposite: the prostate keeps responding, and that response generates diagnoses rather than demonstrated disease. The clinical conclusion converges while the mechanism differs. A modest PSA rise on therapy is expected pharmacodynamics and should be interpreted against a treatment specific baseline rather than population thresholds.

The synthesis is conditional by construction. A reader who rejects the mediation estimate reaches the face value posterior of 1.14 rather than 0.99, and a 12% probability of meaningful harm rather than 1%. The strength of the detection argument is therefore inseparable from the validity of conditioning on PSA, which the negative control does not establish. Two findings do not depend on it: the PSA response to testosterone across the sampled range, and the trial observation that therapy doubles biopsies without increasing cancer.

What we cannot say is that testosterone is inert for the prostate. The organ responds across the physiological range, and in men who already have androgen dependent cancer that response matters, which is why androgen deprivation remains effective therapy.

**Limitations.** Table 6 states each limitation with the expected direction of bias; the full table is in Supplementary Table S9. The most consequential are that lethal disease is untested and would need about 40,000 fatal cases for an equivalence test; that the mediated proportion is not identified, so the channel share is a range and an upper bound; that the collider alternative is not excluded; that the contributing trials excluded men with PSA above about 3.0 and men at elevated prostate cancer risk, which is precisely the stratum in which a PSA triggered detection channel would most plausibly operate, so trial based reassurance does not transport to them; and that MR RAPS is unstable on these data irrespective of implementation, returning 0.009, 0.088 and 0.176 under L2, Huber and Tukey losses (Supplementary Table S6), which is expected of a robust profile estimator at an I squared of 90% and means no single setting should be read as the MR RAPS estimate.

**Conclusion**

Testosterone raises PSA, and much of its apparent cancer risk runs through that channel. For the question clinicians face, whether raising testosterone in men without prostate cancer causes clinically meaningful disease, this study does not establish an answer of no. What it establishes is that a substantial part of the genetic signal runs through PSA, that the residual direct effect is not distinguishable from the null, and that randomized trials over about two years show no excess cancer while doubling biopsies. An absence claim would require an equivalence test that available data cannot support. Prostate safety arguments should rest on randomized outcomes and on tissue sensitivity rather than on receptor saturation, within the population the trials enrolled. The open question is lethal disease, and it remains the priority for future work.

**Declarations**

**Ethical approval.** This study used fully deidentified, publicly available summary statistics from completed genetic association studies and published trial reports. It did not involve human subjects, identifiable private information, or any intervention. The author's institution has been asked to confirm the applicable determination; pending that confirmation the study is described as not constituting human subjects research rather than as carrying an institutional exemption determination.

**Author contributions.** Julian Yin Vieira Borges (Conceptualization, Methodology, Software, Formal analysis, Data curation, Visualization, Writing original draft, Writing review and editing).

**Supplementary data.** Supplementary material, comprising eight Methods sections and Tables S1 to S10, is available at the journal's website.

**Conflict of interest.** None declared.

**Funding.** No specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Data availability.** All data analysed were obtained from publicly accessible repositories: the GWAS Catalog (GCST90012103, GCST90012109, GCST90012113, GCST90461907, GCST006085, GCST002783, GCST90019511, GCST90002238), FinnGen release 12, and the MRC IEU 1000 Genomes European reference panel. Randomized trial counts were taken from published articles and their open supplementary materials. Code, committed result tables, harmonised variant level inputs and figure data are at github.com/julian-borges-md/testosterone-psa-artifact under an MIT licence, at tag v6.1, with that tag archived at https://doi.org/10.5281/zenodo.23324744. No controlled access or patient identifiable data were used or redistributed.

**References**

1. Huggins C, Stevens RE Jr, Hodges CV. Studies on prostatic cancer: II. The effects of castration on advanced carcinoma of the prostate gland. Arch Surg 1941;43:209-23.
2. Morgentaler A. Testosterone and prostate cancer: an historical perspective on a modern myth. Eur Urol 2006;50:935-9.
3. Endogenous Hormones and Prostate Cancer Collaborative Group. Endogenous sex hormones and prostate cancer: a collaborative analysis of 18 prospective studies. J Natl Cancer Inst 2008;100:170-83.
4. Watts EL, Appleby PN, Perez-Cornago A, et al. Low free testosterone and prostate cancer risk: a collaborative analysis of 20 prospective studies. Eur Urol 2018;74:585-94.
5. Bhasin S, Travison TG, Pencina KM, et al. Prostate safety events during testosterone replacement therapy in men with hypogonadism: a randomized clinical trial. JAMA Netw Open 2023;6:e2348692.
6. Lincoff AM, Bhasin S, Flevaris P, et al. Cardiovascular safety of testosterone-replacement therapy. N Engl J Med 2023;389:107-17.
7. Rhomberg-Kauert AM, Tulstrup M, Badsted C, Horwitz H, Heerfordt IM. Mendelian randomization studies of testosterone exposure: a systematic review. Steroids 2025;218:109615.
8. Morgentaler A, Traish AM. Shifting the paradigm of testosterone and prostate cancer: the saturation model and the limits of androgen-dependent growth. Eur Urol 2009;55:310-20.
9. Bañez LL, Hamilton RJ, Partin AW, et al. Obesity-related plasma hemodilution and PSA concentration among men with prostate cancer. JAMA 2007;298:2275-80.
10. Ruth KS, Day FR, Tyrrell J, et al. Using human genetics to understand the disease impacts of testosterone in men and women. Nat Med 2020;26:252-8.
11. Hoffmann TJ, Graff RE, Madduri RK, et al. Genome-wide association study of prostate-specific antigen levels in 392,522 men identifies new loci and improves prediction across ancestry groups. Nat Genet 2025;57:334-44.
12. Schumacher FR, Al Olama AA, Berndt SI, et al. Association analyses of more than 140,000 men identify 63 new prostate cancer susceptibility loci. Nat Genet 2018;50:928-36.
13. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature 2023;613:508-18.
14. Locke AE, Kahali B, Berndt SI, et al. Genetic studies of body mass index yield new insights for obesity biology. Nature 2015;518:197-206.
15. Sinnott-Armstrong N, Tanigawa Y, Amar D, et al. Genetics of 35 blood and urine biomarkers in the UK Biobank. Nat Genet 2021;53:185-94.
16. Chen J, Spracklen CN, Marenne G, et al. The trans-ancestral genomic architecture of glycemic traits. Nat Genet 2021;53:840-60.
17. García-Becerra CA, Arias-Gallardo MI, Juárez-García JE, et al. Cardiovascular and prostate cancer risk associated to testosterone replacement therapy: a systematic review and meta-analysis of 41 randomized controlled trials. Int J Impot Res 2026.
18. Sharifan A. Comment on: cardiovascular and prostate cancer risk associated to testosterone replacement therapy. Int J Impot Res 2026.
19. Zeng B, Qiu S, Xiong X, et al. The effect of different administrations of testosterone therapy on adverse prostate events: a Bayesian network meta-analysis. Front Endocrinol 2022;13:1009900.
20. Burgess S, Butterworth A, Thompson SG. Mendelian randomization analysis with multiple genetic variants using summarized data. Genet Epidemiol 2013;37:658-65.
21. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. Int J Epidemiol 2015;44:512-25.
22. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator. Genet Epidemiol 2016;40:304-14.
23. Hartwig FP, Davey Smith G, Bowden J. Robust inference in summary data Mendelian randomization via the zero modal pleiotropy assumption. Int J Epidemiol 2017;46:1985-98.
24. Xue H, Shen X, Pan W. Constrained maximum likelihood-based Mendelian randomization robust to both correlated and uncorrelated pleiotropic effects. Am J Hum Genet 2021;108:1251-69.
25. Zhao Q, Wang J, Hemani G, Bowden J, Small DS. Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score. Ann Stat 2020;48:1742-69.
26. Morrison J, Knoblauch N, Marcus JH, Stephens M, He X. Mendelian randomization accounting for correlated and uncorrelated pleiotropic effects using genome-wide summary statistics. Nat Genet 2020;52:740-7.
27. Hemani G, Tilling K, Davey Smith G. Orienting the causal relationship between imprecisely measured traits using GWAS summary data. PLoS Genet 2017;13:e1007081.
28. Sanderson E, Spiller W, Bowden J. Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization. Stat Med 2021;40:5434-52.
29. Matthews RAJ. Beyond significance: principles and practice of the Analysis of Credibility. R Soc Open Sci 2018;5:171047.
30. Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: the STROBE-MR Statement. JAMA 2021;326:1614-21.
31. Hynes MC, Watling CZ, Dunneram Y, Key TJ, Perez-Cornago A. Associations of body composition measures with circulating insulin-like growth factor-I, testosterone, and sex hormone-binding globulin concentrations in 16,000 men. Int J Obes 2024;48:1809-17.
32. Vermeulen A, Verdonck L, Kaufman JM. A critical evaluation of simple methods for the estimation of free testosterone in serum. J Clin Endocrinol Metab 1999;84:3666-72.
33. Cheung-Flynn J, Prapapanich V, Cox MB, Riggs DL, Suarez-Quian C, Smith DF. Physiological role for the cochaperone FKBP52 in androgen receptor signaling. Mol Endocrinol 2005;19:1654-66.


**Tables**

**Table 1. Data sources.** Ancestry, units and denominators for every input. Exposure GWAS are male specific.

| Role | Trait and units | Source, ancestry | Sample |
|---|---|---|---|
| Exposure | Bioavailable testosterone, per SD | UK Biobank [10]; European | 178,782 |
| Exposure | SHBG, per SD | UK Biobank [10]; European | 180,726 |
| Exposure | Total testosterone, per SD | UK Biobank [10]; European | 194,453 |
| Bioassay, mediator | PSA, natural log ng/mL | Hoffmann [11]; multiancestry, 28.8 percent non European, 96.5 percent Million Veteran Program | 296,754 discovery; 392,522 with replication |
| Outcome | Prostate cancer, log odds | PRACTICAL [12]; European | 79,148 cases, 61,106 controls (case fraction 0.564) |
| Outcome | Prostate cancer, log odds | FinnGen R12 [13]; Finnish, registry ascertained | 20,368 cases, 156,671 controls (0.115) |
| Covariate | BMI, per SD | GIANT [14]; mixed | 233,572 |
| Covariate | IGF-1, units not stated in deposited metadata | UK Biobank [15]; European | about 350,000 |
| Covariate | Fasting insulin, raw trait scale | MAGIC [16]; trans ancestry | 151,013 |
| Trial evidence | Prostate cancer events | 41 trials identified [17], 10 contributing events, counts from the correspondence supplement [18] | 11,161 randomized; 8,631 in the 10 trials, 50 events |
| Trial evidence | Biopsies and cancer | Zeng [19], randomized trial only stratum | 6,361 |

**Table 2. Complete estimator matrix: effect on PSA, natural log units per SD of exposure.** Bold indicates a 95 percent interval excluding the null. No estimator computed for any exposure is omitted, and all seven are tested in the contrast of Table 3. The MR RAPS row is from the published implementation, mr.raps 0.2, following the cross validation in Supplementary Table S5.

| Estimator | Bioavailable testosterone | Total testosterone | SHBG |
|---|---|---|---|
| IVW | **0.085 (0.042 to 0.127)** | 0.010 (-0.019 to 0.039) | **-0.076 (-0.131 to -0.021)** |
| MR Egger | -0.036 (-0.109 to 0.038) | **-0.071 (-0.120 to -0.021)** | **-0.132 (-0.232 to -0.032)** |
| Weighted median | **0.130 (0.101 to 0.159)** | 0.011 (-0.014 to 0.037) | **-0.069 (-0.111 to -0.028)** |
| Weighted mode | **0.149 (0.102 to 0.196)** | 0.010 (-0.063 to 0.083) | -0.081 (-0.206 to 0.045) |
| IVW, Q outliers removed | **0.145 (0.121 to 0.168)** | **0.038 (0.016 to 0.060)** | **-0.067 (-0.103 to -0.031)** |
| MR cML BIC, data perturbation SE | **0.178 (0.151 to 0.205)** | **0.057 (0.009 to 0.105)** | **-0.060 (-0.113 to -0.007)** |
| MR RAPS Huber (mr.raps 0.2) | **0.088 (0.030 to 0.145)** | -0.020 (-0.068 to 0.028) | **-0.071 (-0.135 to -0.006)** |
| Instruments analysed | 87 | 173 | 232 |
| Cochran Q (I squared) | 859 (90.0%) | 1621 (89.4%) | 1906 (87.9%) |
| MR RAPS tau squared | 2.00e-04 | 2.38e-04 | 8.64e-05 |

**Table 3. Contrast between bioavailable and total testosterone on PSA, all seven estimators.** The difference is tested treating the two estimates as independent, which is conservative because they share the PSA sample. An asterisk marks the row where both estimates are negative, so the comparison is not a magnitude contrast between positive effects.

| Estimator | Bioavailable | Total | Difference (SE) | z | p |
|---|---|---|---|---|---|
| IVW | 0.0848 | 0.0104 | +0.0743 (0.0263) | +2.83 | 0.0047 |
| MR Egger * | -0.0359 | -0.0707 | +0.0349 (0.0451) | +0.77 | 0.4401 |
| Weighted median | 0.1298 | 0.0113 | +0.1185 (0.0197) | +6.03 | <0.0001 |
| Weighted mode | 0.1490 | 0.0096 | +0.1394 (0.0443) | +3.14 | 0.0017 |
| IVW, Q outliers removed | 0.1447 | 0.0380 | +0.1067 (0.0162) | +6.59 | <0.0001 |
| MR cML BIC, data perturbation SE | 0.1780 | 0.0569 | +0.1211 (0.0280) | +4.33 | <0.0001 |
| MR RAPS Huber (mr.raps 0.2) | 0.0876 | -0.0200 | +0.1076 (0.0381) | +2.83 | 0.0047 |

**Table 4. Complete estimator matrix: prostate cancer, odds ratio per SD of exposure.** Bold indicates a 95 percent interval excluding unity. No estimator computed for any exposure and cohort is omitted.

| Estimator | Bioavailable T, PRACTICAL | Total T, PRACTICAL | SHBG, PRACTICAL | Bioavailable T, FinnGen |
|---|---|---|---|---|
| IVW | **1.193 (1.079 to 1.318)** | 0.998 (0.913 to 1.092) | **0.838 (0.709 to 0.991)** | 1.054 (0.906 to 1.226) |
| MR Egger | 1.097 (0.906 to 1.328) | 1.020 (0.870 to 1.195) | 0.945 (0.700 to 1.275) | 0.845 (0.618 to 1.156) |
| Weighted median | **1.127 (1.011 to 1.256)** | 1.057 (0.980 to 1.141) | 0.884 (0.760 to 1.027) | 0.984 (0.842 to 1.149) |
| Weighted mode | 1.108 (0.988 to 1.242) | 1.063 (0.981 to 1.151) | 0.874 (0.707 to 1.079) | 0.960 (0.801 to 1.150) |
| IVW, Q outliers removed | **1.153 (1.072 to 1.239)** | **1.075 (1.011 to 1.144)** | 0.922 (0.825 to 1.032) | 0.974 (0.855 to 1.110) |
| MR cML BIC, data perturbation SE | **1.155 (1.062 to 1.256)** | **1.088 (1.019 to 1.161)** | 0.956 (0.830 to 1.102) | 0.961 (0.844 to 1.095) |
| MR RAPS Huber (mr.raps 0.2) | **1.160 (1.056 to 1.275)** | 1.064 (0.977 to 1.158) | 0.979 (0.843 to 1.137) | 0.969 (0.815 to 1.150) |
| Instruments analysed | 87 | 174 | 232 | 82 |
| Cochran Q (I squared) | 232 (62.9%) | 790 (78.1%) | 940 (75.4%) | 213 (62.0%) |
| MR RAPS tau squared | 1.84e-04 | 4.20e-04 | 2.68e-04 | 6.69e-04 |

**Table 5. Posterior odds ratio for any prostate cancer.** The likelihood is the exact conditional likelihood over the 10 trials reporting prostate cancer events (8,631 participants, 50 events), not over all 41 trials identified by the source review.

| Prior | Posterior OR (95 percent interval) | P(OR above 1.2) |
|---|---|---|
| None (vague) | 0.89 (0.50 to 1.56) | 15 percent |
| Genetics at face value, pooled cohorts | 1.14 (1.05 to 1.24) | 12 percent |
| Genetics, PSA channel removed, pooled cohorts | 0.99 (0.84 to 1.17) | 1 percent |
| Genetics at face value, PRACTICAL only | 1.18 (1.07 to 1.30) | 38 percent |
| Genetics, PSA channel removed, PRACTICAL only | 1.02 (0.83 to 1.25) | 6 percent |
| Genetics, PSA channel removed, pooled, largest trial excluded | 0.99 (0.83 to 1.17) | 1 percent |

**Table 6. Principal limitations and the expected direction of bias.** The full 23 row table is Supplementary Table S9.

| Issue | Consequence and expected direction |
|---|---|
| Lethal disease untested | The clinically decisive question remains open; equivalence within OR 1.10 would need about 40,000 fatal cases, at which point the two one sided test reaches 81 percent power |
| The mediated proportion is not identified | The odds ratio is not collapsible at a case fraction of 0.564 and the additivity gap is 27.3 percent of the total effect, so the channel share is a product method range and an upper bound, not a point estimate |
| The collider alternative is not excluded | Two of three negative controls are uninformative in the available instrument sets, so a collider component cannot be ruled out and the detection share may be overstated |
| Trial evidence does not transport to men near a biopsy threshold | The contributing trials excluded men with PSA above about 3.0 and men at elevated prostate cancer risk, which is the stratum in which the detection channel would most plausibly operate |
| A linear IV average cannot falsify a threshold model | A positive average is compatible with saturation if enough men lie below threshold; we bound this but do not test threshold location, and the bound assumes population share weighting while the estimator uses inverse variance weighting |
| MR RAPS is unstable on these data | The published implementation returns 0.009, 0.088 and 0.176 for the primary analysis under L2, Huber and Tukey losses with overdispersion enabled, so no single setting should be read as the MR RAPS estimate |
| Directional pleiotropy and heterogeneity | Egger intercept p = 1.3e-4 and I squared 87.9 to 90.0 percent in the PSA analyses, traced in part to tissue sensitivity loci, biasing IVW toward the null |
| Mixed ancestry across the analysis | European exposures against a multiancestry PSA outcome, a European cancer cohort and a Finnish one; no transportability claim is made to African ancestry men, who carry the highest disease burden |
| Prespecification is internal | Time stamped in the committed analysis log, not registered with an external registry |

**Figure legends**

**Figure 1. The complete estimator matrix for PSA.** Effect of each exposure on PSA in natural log units per standard deviation, for all seven estimators and all three exposures. Filled markers indicate a 95 percent interval excluding the null and open markers an interval containing it; bars are 95 percent confidence intervals. Values are taken from the same committed matrix as Table 2.

**Figure 2. Randomized trial evidence and the Bayesian synthesis.** A, testosterone therapy in randomized trials: prostate biopsies and prostate cancer. B, posterior odds ratio for any prostate cancer by assumed therapeutic shift, with and without the PSA channel.

![](Figure1_saturation_full_matrix.png)

![](Figure4_trials_and_synthesis.png)
