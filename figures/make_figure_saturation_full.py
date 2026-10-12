"""Figure 1: the complete estimator matrix for PSA, all seven
estimators for all three exposures, so the figure matches the completeness
commitment made in the Table 2 caption."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300})
C = {'biot': '#1f4e79', 'tt': '#7f7f7f', 'shbg': '#c55a11'}
M = pd.read_csv('results/stage7_matrix/estimator_matrix.tsv', sep='\t')
EST = ['IVW', 'MR Egger', 'Weighted median', 'Weighted mode',
       'IVW, Q outliers removed', 'MR cML BIC (DP SE)', 'MR RAPS Huber']
LAB = {'IVW': 'IVW', 'MR Egger': 'MR Egger', 'Weighted median': 'Weighted median',
       'Weighted mode': 'Weighted mode', 'IVW, Q outliers removed': 'IVW, Q outliers removed',
       'MR cML BIC (DP SE)': 'MR cML (data perturbation SE)', 'MR RAPS Huber': 'MR RAPS (mr.raps 0.2)'}
PANELS = [('biot->psa', 'Bioavailable testosterone', 'biot'),
          ('tt->psa', 'Total testosterone', 'tt'),
          ('shbg->psa', 'SHBG', 'shbg')]

fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.6), sharex=True)
for ax, (an, title, ck) in zip(axes, PANELS):
    sub = M[M.analysis == an].set_index('estimator')
    y = np.arange(len(EST))[::-1]
    for yi, e in zip(y, EST):
        r = sub.loc[e]
        sig = bool(r.excludes_null)
        ax.errorbar(r.beta, yi, xerr=[[r.beta - r.lo], [r.hi - r.beta]], fmt='o',
                    color=C[ck], ms=4.5, capsize=2, lw=1.2,
                    mfc=C[ck] if sig else 'white', mew=1.2)
    ax.axvline(0, color='black', lw=0.8, ls=':')
    ax.set_yticks(y)
    ax.set_yticklabels([LAB[e] for e in EST] if ck == 'biot' else [''] * len(EST))
    ax.set_title(title, fontsize=9)
    ax.set_xlim(-0.26, 0.26)
axes[1].set_xlabel('Effect on PSA (natural log units per SD of exposure, 95% CI)')
axes[1].annotate('Saturation predicts zero', xy=(0, 0.15), xytext=(0.075, 0.55),
                 fontsize=7, style='italic', va='center',
                 arrowprops=dict(arrowstyle='-', lw=0.6, color='0.35'))
plt.tight_layout()
plt.savefig('figures/Figure1_saturation_full_matrix.png', bbox_inches='tight')
plt.savefig('figures/Figure1_saturation_full_matrix.pdf', bbox_inches='tight')
print('wrote figures/Figure1_saturation_full_matrix.{png,pdf}')

rows = []
for an, title, _ in PANELS:
    sub = M[M.analysis == an].set_index('estimator')
    for e in EST:
        r = sub.loc[e]
        rows.append(dict(panel=title, estimator=LAB[e], beta=r.beta, lo=r.lo, hi=r.hi,
                         p=r.p, nsnp=int(r.nsnp), excludes_null=bool(r.excludes_null),
                         source='stage7_matrix/estimator_matrix.tsv'))
pd.DataFrame(rows).to_csv('figures/data/Figure1_saturation_full_matrix_data.tsv',
                          sep='\t', index=False)
print('wrote plotted values, %d rows' % len(rows))
