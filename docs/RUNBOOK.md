# RLO_TPC | Runbook: full reanalysis from raw summary statistics

`make verify` needs none of this. Use this runbook to regenerate the committed results from source.

## 1. Retrieve and verify inputs

```bash
bash data/download_sources.sh          # about 3.3 GB into work/, SHA256 verified
```

Place PLINK 1.90b7.2 and the MRC IEU EUR reference panel in `work/tools/`.

## 2. Point the scripts at the working tree

Every stage script resolves two roots from the environment and falls back to the
paths the original run used, so either route works. Setting the variables needs
no symlink and no write access to /tmp:

```bash
export GW_WORK="$PWD/work"
export GW_TOOLS="$PWD/work/tools"
```

Or keep the defaults and link them, which is what the original run did:

```bash
ln -s "$PWD/work" /tmp/gw
ln -s "$PWD/work/tools" /tmp/tools
```

Then place the scripts where each stage runs:

```bash
cp src/stage1_saturation/*.py work/
cp src/stage2_detection/*.py work/p2/ && cp src/stage2_detection/robust2.py src/stage2_detection/cml_dp.py work/
mkdir -p work/p3 && cp src/stage3_synthesis/*.py results/stage3/trial_counts_from_sharifan_supp.tsv work/p3/
```

## 3. Stage 1: saturation bioassay (work/)

| Step | Command |
|---|---|
| Genome wide significant variants | `python3 sig.py <file> <out> p_value` for each exposure |
| QC and clumping input | `python3 prep.py`, then PLINK `--clump` at r2 0.001, 10 Mb |
| Extract instruments | `python3 extract.py <job>` for biot, shbg, tt, psa, prac, fg |
| Univariable MR | `python3 run1.py` |
| Outlier inspection | `python3 run2.py` |
| Multivariable MR | `python3 run3.py` |
| Colocalization | `python3 region.py` for each exposure, then `python3 coloc.py` |

## 4. Stage 2: rival explanations (work/p2/)

| Step | Command |
|---|---|
| Reverse MR | `python3 extract2.py`, then `python3 run_reverse.py` |
| Metabolic instruments | `python3 sig3.py <bmi or igf or fi>`, then PLINK clumping |
| Metabolic MVMR | `python3 ext_mv.py <file> <out> <mode>`, then `python3 run_mvmr2.py` |
| Robust estimators | from work/: `python3 robust2.py`, then `python3 cml_dp.py` |

## 5. Stage 3: synthesis (work/p3/)

```bash
python3 exact.py
```

## Note on the path change

The scripts were written against `/tmp/gw` and `/tmp/tools`, which was recorded
as limitation L1: a reader had to recreate a scratch location before anything
would run. On 19 October 2026 the two roots became configurable, with those
values kept as the defaults. The substitution was mechanical across 23 files and
84 call sites and changes no analysis, but the committed results were produced
by the code as it stood on 11 October, before the edit. L1 is retired.

## 6. Compare

Copy outputs into `results/` and run `make verify`. Any disagreement prints DRIFT.
