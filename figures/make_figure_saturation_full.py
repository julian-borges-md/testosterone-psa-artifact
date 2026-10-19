"""Figure 1: the complete estimator matrix for PSA, all seven estimators for
all three exposures, so the figure matches the completeness commitment made in
the Table 2 caption.

Laid out as three stacked panels rather than three side by side. At the
Society's two column width of 5.8 inches, three panels in a row leave about
1.5 inches of plot area each, which is not enough for the estimator labels at
8 point. Stacking spends the height allowance instead, which is 9.25 inches
and was previously unused."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style
_style.apply()
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

M = pd.read_csv('results/stage7_matrix/estimator_matrix.tsv', sep='\t')
EST = ['IVW', 'MR Egger', 'Weighted median', 'Weighted mode',
       'IVW, Q outliers removed', 'MR cML BIC (DP SE)', 'MR RAPS Huber']
LAB = {'IVW': 'IVW', 'MR Egger': 'MR Egger', 'Weighted median': 'Weighted median',
       'Weighted mode': 'Weighted mode', 'IVW, Q outliers removed': 'IVW, Q outliers removed',
       'MR cML BIC (DP SE)': 'MR cML (data perturbation SE)', 'MR RAPS Huber': 'MR RAPS (mr.raps 0.2)'}
PANELS = [('biot->psa', 'Bioavailable testosterone', 'biot', 'a'),
          ('tt->psa', 'Total testosterone', 'tt', 'b'),
          ('shbg->psa', 'SHBG', 'shbg', 'c')]

fig, axes = plt.subplots(3, 1, figsize=(_style.TWO_COLUMN, 6.4), sharex=True)
for ax, (an, title, ck, letter) in zip(axes, PANELS):
    sub = M[M.analysis == an].set_index('estimator')
    y = np.arange(len(EST))[::-1]
    for yi, e in zip(y, EST):
        r = sub.loc[e]
        sig = bool(r.excludes_null)
        ax.errorbar(r.beta, yi, xerr=[[r.beta - r.lo], [r.hi - r.beta]], fmt='o',
                    color=_style.C[ck], ms=4.0, capsize=2, lw=1.0,
                    mfc=_style.C[ck] if sig else 'white', mew=1.0)
    ax.axvline(0, color='black', lw=0.8, ls=':')
    ax.set_yticks(y)
    ax.set_yticklabels([LAB[e] for e in EST])
    ax.set_ylim(-0.7, len(EST) - 0.3 + (0.6 if ck == 'tt' else 0.0))
    ax.set_title(title, loc='left')
    ax.set_xlim(-0.26, 0.26)
    _style.panel_label(ax, letter)
axes[2].set_xlabel('Effect on PSA (natural log units per SD of exposure, 95% CI)')
# Placed in the total testosterone panel, which is where the prediction is met,
# and beside the null line rather than on a leader that crosses an interval.
axes[1].text(0.012, len(EST) - 0.5, 'Saturation predicts zero',
             style='italic', va='center', ha='left')
plt.tight_layout()
_style.check(fig, 'Figure1_saturation_full_matrix')
plt.savefig('figures/Figure1_saturation_full_matrix.png', bbox_inches='tight')
plt.savefig('figures/Figure1_saturation_full_matrix.pdf', bbox_inches='tight')
plt.savefig('figures/Figure1_saturation_full_matrix.eps', bbox_inches='tight')
print('wrote figures/Figure1_saturation_full_matrix.{png,pdf,eps}')

rows = []
for an, title, _ck, _l in PANELS:
    sub = M[M.analysis == an].set_index('estimator')
    for e in EST:
        r = sub.loc[e]
        rows.append(dict(panel=title, estimator=LAB[e], beta=r.beta, lo=r.lo, hi=r.hi,
                         p=r.p, nsnp=int(r.nsnp), excludes_null=bool(r.excludes_null),
                         source='stage7_matrix/estimator_matrix.tsv'))
pd.DataFrame(rows).to_csv('figures/data/Figure1_saturation_full_matrix_data.tsv',
                          sep='\t', index=False)
print('wrote plotted values, %d rows' % len(rows))
