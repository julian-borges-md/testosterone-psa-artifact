# Figure alt text and legends

The journal requires alt text for every figure in the main manuscript file,
directly under the figure legend and prefixed `Alt text:`. The manuscript is
therefore the source: these are the same paragraphs, lifted from it, for a
submission system that also asks for alt text in a field of its own.

They are not maintained separately. `src/check_manuscript.py` asserts that each
paragraph here is byte for byte the paragraph in the manuscript, and asserts the
two claims a figure rebuild actually invalidates, the panel layout and the axis
range, against the figure script that writes the file. The alt text written for
the previous venue survived the 19 October figure rebuild describing a layout the
figures no longer had, because nothing read it.

## Figure 1. The complete estimator matrix for PSA

**Alt text.** Three forest plots stacked vertically and sharing one horizontal axis, labelled (a) bioavailable testosterone in dark blue, (b) total testosterone in grey and (c) SHBG in orange. The axis is the effect on PSA in natural log units per standard deviation of exposure, running from minus 0.26 to 0.26 with ticks at minus 0.2, minus 0.1, zero, 0.1 and 0.2 and a dotted vertical line at zero. Each panel lists the same seven estimators down the vertical axis in the same order, from inverse variance weighted at the top through MR Egger, weighted median, weighted mode, inverse variance weighted with Q outliers removed and MR cML to MR RAPS at the bottom. Filled markers mark intervals that exclude zero and open markers intervals that contain it. In panel (a) six of the seven estimates lie to the right of zero with intervals that exclude it, MR Egger alone lying to the left with an interval that crosses zero. Panel (b) carries the annotation Saturation predicts zero, and its estimates cluster near zero: two exclude zero on the positive side, MR Egger excludes it on the negative side, and the remaining four intervals contain it. In panel (c) six of the seven estimates lie to the left of zero with intervals that exclude it, the weighted mode alone crossing it. Values are listed in Table 2.

**Legend as printed.** Effect of each exposure on PSA in natural log units per standard deviation, for all seven estimators and all three exposures. Filled markers indicate a 95 percent interval excluding the null and open markers an interval containing it; bars are 95 percent confidence intervals. Values are taken from the same committed matrix as Table 2.

## Figure 2. Randomized trial evidence and the Bayesian synthesis

**Alt text.** Two panels side by side. Panel (a) is a forest plot of relative risks from randomized trials of testosterone therapy, on a logarithmic horizontal axis labelled 0.5, 1, 2 and 5 with a dotted vertical line at 1. Prostate biopsies, an orange square, sit above 1 at about 2.4 with an interval running from just above 1 to beyond 5. Prostate cancer, a dark blue square, sits below 1 at about 0.85 with an interval from about 0.44 to about 1.6 that spans 1. Panel (b) plots the posterior odds ratio for any prostate cancer, on a vertical axis from 0.6 to 1.8, against the assumed therapeutic shift in standard deviations, from 0.5 to 3.0. A dotted horizontal line marks 1 and a black line of alternating dashes and dots marks 1.2, labelled meaningful harm (OR 1.2). A solid orange line with circular markers, for the genetic evidence taken at face value, rises from about 1.07 to about 1.39 and crosses the 1.2 line at a shift of about 1.4. A dashed dark blue line with square markers, for the genetic evidence with the PSA channel removed, stays flat and slightly below 1 across the whole range. Both curves carry shaded 95 percent credible bands that widen as the assumed shift grows; the band around the dashed blue curve spans values well below and well above 1 throughout.

**Legend as printed.** (a), testosterone therapy in randomized trials: prostate biopsies and prostate cancer. (b), posterior odds ratio for any prostate cancer by assumed therapeutic shift, with and without the PSA channel.
