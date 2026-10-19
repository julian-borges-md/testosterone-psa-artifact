# Cover letter

Julian Yin Vieira Borges, MD, MS
Department of Computer Science, Boston University Metropolitan College
Boston, MA 02215, United States
jyborges@bu.edu | ORCID 0009-0001-9929-3135

12 October 2026

Dear Editors,

I am submitting "Testosterone Raises PSA Across the Physiological Range, with Implications for Detection Bias in Prostate Cancer" for consideration as a Research Article.

The mechanism is the premise of the design rather than an inference from it. Prostate specific antigen is transcribed under direct androgen receptor control, which is what makes it an in vivo bioassay of prostatic androgen signalling rather than a correlate of one, and the completed STROBE-MR checklist accompanies the submission with every item addressed.

For decades clinical practice has carried the belief that testosterone causes prostate cancer. The usual modern reconciliation of the conflicting evidence is the androgen saturation model, which holds that prostatic androgen receptors are saturated at low serum testosterone, so variation within the physiological range cannot affect the prostate. That model is widely cited in clinical guidance and it has never been tested causally in humans.

This study tests the prediction it makes. Because prostate specific antigen is transcribed under direct androgen receptor control, PSA serves as an in vivo bioassay of prostatic androgen signalling: if receptors were saturated, genetic variation in bioavailable testosterone within the physiological range should not move PSA. It does. Using male specific instruments in 178,782 to 194,453 men and PSA in 392,522, genetically higher bioavailable testosterone raised PSA by 0.085 natural log units per standard deviation, with intervals excluding the null in six of seven estimators.

I want to be precise about what that does and does not establish, because an earlier version of this manuscript was not. A linear instrumental variable average over a population in which some men lie below the threshold is compatible with strict saturation, so this design bounds rather than locates a threshold and cannot falsify the model. The bound is informative: for the observed average to arise solely from men below threshold, their slope would have to be about 11.5 times the average at 8 nmol/L and about 15,000 times at 4 nmol/L. It is still a bound, and the manuscript says so in its abstract, its Results and its limitations table.

The second finding is the one in the title. The prostate keeps responding to testosterone; what that response produces is PSA, and a raised PSA brings forward the biopsy that finds a tumour already present. Between a third and a half of the genetic association with diagnosed prostate cancer is carried by the PSA channel, and conditioning on PSA reduces the direct effect to an odds ratio of 1.04 with an interval from 0.83 to 1.30. Randomized trials show the same channel directly: therapy more than doubled prostate biopsies without increasing cancer. Reverse causation and confounding by BMI, IGF-1 and fasting insulin were excluded.

Combining both bodies of evidence in a Bayesian synthesis, remaining 95 percent convinced that testosterone therapy causes prostate cancer requires beginning from a belief that risk is at least doubled, held as firmly as the randomized evidence itself. With the PSA channel removed the posterior odds ratio is 0.99.

The readership this most concerns is the one that prescribes testosterone and then reads a PSA on therapy. A modest rise is expected pharmacodynamics on this evidence rather than a signal of tumour biology, and it should be interpreted against a treatment specific baseline. That is a statement about how to read a biomarker, not a safety clearance, and the manuscript is careful to keep the two apart.

Three features may interest your reviewers.

First, the reproducibility apparatus is unusually complete, and it is adversarial rather than decorative. Hypotheses, priors and decision rules were time stamped before estimation. The public repository carries the harmonised variant level inputs, so every cell of the estimator matrices is recomputable offline. Five assertion modules re-derive every reported number and exit non-zero on any disagreement; one of them exists only to prove that another can fail. Every estimator has been cross validated against a published package on identical inputs, and where that comparison disagreed with this study's own code the published implementation was adopted and the change was reported, including the two counts it moved in the manuscript's favour.

Second, the study reaches the clinical conclusion favoured by proponents of testosterone therapy through the opposite mechanism from the one they invoke, which makes the argument harder to attribute to prior commitment.

Third, the limitations are stated with the expected direction of bias for every methodological choice, including the two that most constrain the conclusion. Conditioning on PSA may carry a collider component, so the channel share is reported as a range with an upper bound rather than as a point estimate. And MR RAPS is unstable on this instrument set irrespective of implementation, so every count that rests on that row states the loss setting it depends on and the alternative under a different setting.

What the study does not show is stated as plainly. No public data set has enough fatal cases to test whether androgens affect lethal prostate cancer; detecting an odds ratio of 1.19 for aggressive disease would need about 40,000 fatal cases. That question is named as the priority for future work rather than answered here.

The manuscript is original, is not under consideration elsewhere and has not been published previously. All data analysed are publicly available and no controlled access or patient identifiable data were used or redistributed. The analysis code is openly released under an MIT licence and archived at https://doi.org/10.5281/zenodo.23483135. I have no conflicts of interest to declare and received no funding for this work. The completed STROBE-MR checklist accompanies the submission.

Thank you for considering this manuscript.

Yours sincerely,

Julian Yin Vieira Borges, MD, MS
