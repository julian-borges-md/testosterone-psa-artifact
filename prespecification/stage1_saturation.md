# RLO_TPC_SaturationMR_v1.0 prespecified decision rules
Logged: 2026-10-09T17:13:24Z, before any exposure to outcome estimate was computed.

Primary test (saturation bioassay): MR of bioavailable testosterone (men, per SD) on PSA.
- Saturation SUPPORTED if the IVW 95% CI includes zero AND no robust estimator (weighted median, weighted mode, Egger) shows a positive estimate with 95% CI excluding zero.
- Saturation FALSIFIED if the IVW estimate is positive with 95% CI excluding zero AND at least two of the three robust estimators agree in direction with CI excluding zero, after removal of Q outliers.
- Otherwise INCONCLUSIVE.
Free hormone contrast: if SHBG raises total testosterone but not bioavailable testosterone, free hormone hypothesis predicts SHBG to PSA effect near zero or negative.
Instrument rules: p<5e-8, r2<0.001 in 10 Mb, 1000G EUR (n=503), biallelic SNPs, MAF>=0.01, INFO>=0.8 (UKB exposures), palindromic SNPs with MAF>0.42 removed, EAF discordance >0.2 removed.
Steiger: observed scale (case fraction in PRACTICAL >> population prevalence; lab rule).
