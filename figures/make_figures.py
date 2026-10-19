import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style
_style.apply()
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import NullFormatter, NullLocator
C = _style.C
R=pd.read_csv('results/stage1/mr_results.tsv',sep='\t'); RB=pd.read_csv('results/stage2/robust_cml_raps.tsv',sep='\t'); DP=pd.read_csv('results/stage2/cml_dp.tsv',sep='\t')
U=pd.read_csv('results/stage2/uni_metabolic.tsv',sep='\t'); UP=pd.read_csv('results/stage2/uni_metabolic_psa.tsv',sep='\t')
S=pd.read_csv('results/stage4/synthesis_pooled.tsv',sep='\t')
MX=pd.read_csv('results/stage7_matrix/estimator_matrix.tsv',sep='\t')
def mx(a,e):
    r=MX[(MX.analysis==a)&(MX.estimator==e)].iloc[0]; return r.beta,r.lo,r.hi
def get(df,a,m): r=df[(df.analysis==a)&(df.method==m)].iloc[0]; return r.beta,r.lo,r.hi

# ---------- Figure 1: conceptual model ----------
fig,ax=plt.subplots(figsize=(_style.TWO_COLUMN,2.9)); ax.axis('off'); ax.set_xlim(0,10.4); ax.set_ylim(0,5)
def box(x,y,t,fc): ax.add_patch(FancyBboxPatch((x-0.95,y-0.38),1.9,0.76,boxstyle='round,pad=0.05',fc=fc,ec='black',lw=0.8)); ax.text(x,y,t,ha='center',va='center')
def arr(a,b,t='',c='black',ls='-',dy=0.15,rad=None):
    ap=dict(arrowstyle='->',color=c,lw=1.2,ls=ls)
    if rad is not None: ap['connectionstyle']='arc3,rad=%s'%rad
    ax.annotate('',xy=b,xytext=a,arrowprops=ap)
    if t: ax.text((a[0]+b[0])/2,(a[1]+b[1])/2+dy,t,ha='center',color=c)
box(1.1,3.5,'Bioavailable\ntestosterone','#dbe5f1'); box(3.6,3.5,'Prostate AR\nsignalling','#dbe5f1')
box(6.4,4.3,'PSA','#fde9d9'); box(9.2,4.3,'Biopsy','#fde9d9'); box(9.2,2.6,'Diagnosed\nprostate cancer','#f2f2f2')
box(6.4,1.4,'Tumour initiation\nor progression','#e2efda'); box(9.2,1.0,'Lethal prostate\ncancer','#e2efda')
box(1.1,1.6,'BMI\n(hemodilution)','#ebf1de')
arr((2.05,3.5),(2.65,3.5)); arr((4.55,3.75),(5.45,4.2)); ax.text(4.75,4.35,'raises PSA\n(saturation bounded)',ha='center',color=C['biot'])
arr((7.35,4.3),(8.25,4.3),c=_style.PATH_DETECTION); arr((9.2,3.92),(9.2,2.98),c=_style.PATH_DETECTION); ax.text(9.35,3.45,'detection',color=_style.PATH_DETECTION)
arr((4.4,3.12),(5.6,1.75),c=_style.PATH_BIOLOGY,ls='--'); ax.text(4.55,2.3,'biology?',color=_style.PATH_BIOLOGY)
arr((7.35,1.6),(8.25,2.35),c=_style.PATH_BIOLOGY,ls='--'); arr((7.35,1.2),(8.25,1.0),c=_style.PATH_BIOLOGY,ls='--')
arr((2.05,1.6),(5.45,3.95),c=C['tt'],rad=-0.28); ax.text(3.6,2.1,'lowers PSA',color=C['tt'],rotation=22,ha='center')
ax.text(5.2,0.1,'Solid path: PSA driven detection (supported). Dashed path: tumour biology (not demonstrated for lethal disease).',ha='center',style='italic')
_style.check(fig,'Figure1_conceptual_model'); plt.savefig('figures/Figure1_conceptual_model.png',bbox_inches='tight'); plt.savefig('figures/Figure1_conceptual_model.pdf',bbox_inches='tight'); plt.close()

# ---------- Figure 2: saturation and free hormone test ----------
rows=[('Bioavailable T, IVW',*mx('biot->psa','IVW'),C['biot']),
      ('Bioavailable T, weighted median',*mx('biot->psa','Weighted median'),C['biot']),
      ('Bioavailable T, weighted mode',*mx('biot->psa','Weighted mode'),C['biot']),
      ('Bioavailable T, IVW Q outliers removed',*mx('biot->psa','IVW, Q outliers removed'),C['biot']),
      ('Bioavailable T, MR cML (DP SE)',*mx('biot->psa','MR cML BIC (DP SE)'),C['biot']),
      ('Bioavailable T, MR RAPS (mr.raps 0.2)',*mx('biot->psa','MR RAPS Huber'),C['biot']),
      ('Bioavailable T, MR Egger',*mx('biot->psa','MR Egger'),C['biot']),
      ('Total T, IVW',*mx('tt->psa','IVW'),C['tt']),
      ('SHBG, IVW',*mx('shbg->psa','IVW'),C['shbg'])]
fig,ax=plt.subplots(figsize=(_style.TWO_COLUMN,4.0)); y=np.arange(len(rows))[::-1]
for yi,(lab,b,lo,hi,c) in zip(y,rows):
    ax.errorbar(b,yi,xerr=[[b-lo],[hi-b]],fmt='o',color=c,ms=5,capsize=2,lw=1.2)
ax.axvline(0,color='black',lw=0.8,ls=':')
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.set_xlabel('Effect on PSA (natural log units per SD of exposure, 95% CI)')
ax.text(0.003,len(rows)-0.4,'Saturation predicts zero',style='italic')
plt.tight_layout(); _style.check(fig,'Figure2_saturation_test'); plt.savefig('figures/Figure2_saturation_test.png'); plt.savefig('figures/Figure2_saturation_test.pdf'); plt.close()

# ---------- Figure 3: detection signature ----------
ex=[('Bioavailable\ntestosterone','biot',C['biot']),('BMI','bmi',C['bmi']),('IGF-1','igf',C['igf'])]
fig,axs=plt.subplots(1,3,figsize=(_style.TWO_COLUMN,2.1),sharey=True)
titles=['Effect on PSA\n(log units)','Prostate cancer, PRACTICAL\n(OR; heavily screened)','Prostate cancer, FinnGen\n(OR; registry)']
for j,ax in enumerate(axs):
    for i,(lab,k,c) in enumerate(ex):
        if j==0: src=R if k=='biot' else UP; a=f'{k}->psa'
        else:
            o='prac' if j==1 else 'fg'; src=R if k=='biot' else U; a=f'{k}->{o}'
        b,lo,hi=get(src,a,'IVW')
        if j>0: b,lo,hi=np.exp([b,lo,hi])
        yi=2-i; ax.errorbar(b,yi,xerr=[[b-lo],[hi-b]],fmt='o',color=c,ms=5,capsize=2,lw=1.2)
    ax.axvline(0 if j==0 else 1,color='black',lw=0.8,ls=':'); ax.set_title(titles[j],pad=4)
    _style.panel_label(ax,'abc'[j],x=0.0,y=1.30,ha='left')
    ax.set_ylim(-0.55,2.55)
axs[0].set_yticks([2,1,0]); axs[0].set_yticklabels([e[0] for e in ex])
# Panels b and c are both odds ratios and are read against each other, so they
# carry one scale. Different ranges invited a visual comparison the axes did
# not support, which is the comparison the between source analysis now makes
# formally.
_or = min(axs[1].get_xlim()[0], axs[2].get_xlim()[0]), max(axs[1].get_xlim()[1], axs[2].get_xlim()[1])
axs[1].set_xlim(_or); axs[2].set_xlim(_or)
plt.tight_layout(); _style.check(fig,'Figure3_detection_signature'); plt.savefig('figures/Figure3_detection_signature.png'); plt.savefig('figures/Figure3_detection_signature.pdf'); plt.close()

# ---------- Figure 4: trials and synthesis ----------
fig,axs=plt.subplots(1,2,figsize=(_style.TWO_COLUMN,3.1),gridspec_kw={'width_ratios':[1,1.3]})
# Colour. The greyscale here was a print colour charge decision taken for a
# venue with a print edition. The current venue is online only and fully open
# access, so the charge does not exist. Line style and marker shape still carry
# the distinction, so the panel survives a monochrome print by a reader.
ax=axs[0]; tr=[('Prostate biopsies',2.38,1.01,5.60,_style.PATH_DETECTION),('Prostate cancer',0.85,0.44,1.62,_style.C['biot'])]
for yi,(lab,b,lo,hi,c) in zip([1,0],tr): ax.errorbar(b,yi,xerr=[[b-lo],[hi-b]],fmt='s',color=c,ms=6,capsize=3,lw=1.4)
ax.axvline(1,color='black',lw=0.8,ls=':'); ax.set_xscale('log')
# On a log axis spanning about one decade matplotlib also labels the MINOR
# ticks, with LogFormatterSciNotation, and those labels overprint the major
# ones: the published panel read '4x10^-1 6x10^-1 ... 2 3x10^0 4x10^0' on top
# of 0.5, 1, 2 and 5. Both the minor locator and its formatter are therefore
# removed, leaving only the four explicit major labels.
ax.set_xticks([0.5,1,2,5]); ax.set_xticklabels(['0.5','1','2','5'])
ax.xaxis.set_minor_locator(NullLocator()); ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlim(0.38, 6.6)
ax.set_yticks([1,0]); ax.set_yticklabels([t[0] for t in tr]); ax.set_ylim(-0.6,1.6); ax.set_xlabel('Relative risk, randomized trials (95% CI)'); ax.set_title('Testosterone therapy in trials',loc='left'); _style.panel_label(ax,'a')
ax=axs[1]; d=S
for pr,c,sty,lab in [('face value, pooled PRACTICAL+FinnGen',_style.C['shbg'],'-o','Genetics at face value'),('PSA channel removed, pooled',_style.C['biot'],'--s','Genetics, PSA channel removed')]:
    # Solid tints rather than alpha. The PostScript backend does not support
    # transparency and renders a partially transparent band opaque, which in the
    # submitted EPS would have turned these bands into solid blocks covering the
    # curves. The tint is the same colour the alpha produced, computed once.
    q=d[d.prior==pr]
    ax.fill_between(q.shift_SD,q.lo,q.hi,color=_style.tint(c,0.15),lw=0,zorder=1)
    ax.plot(q.shift_SD,q.post_OR,sty,color=c,ms=3,label=lab,zorder=3)
ax.axhline(1,color='black',lw=0.8,ls=':'); ax.axhline(1.2,color='black',lw=0.7,ls='-.'); ax.text(0.55,1.22,'meaningful harm (OR 1.2)',color='black',
        bbox=dict(facecolor='white',edgecolor='none',pad=1.0))
ax.set_xlabel('Assumed therapy shift (SD)'); ax.set_ylabel('Posterior OR, any prostate cancer'); ax.legend(frameon=False,loc='upper left')
ax.set_title('Bayesian synthesis with 41 trials',loc='left'); _style.panel_label(ax,'b')
plt.tight_layout(); _style.check(fig,'Figure4_trials_and_synthesis'); plt.savefig('figures/Figure4_trials_and_synthesis.png'); plt.savefig('figures/Figure4_trials_and_synthesis.pdf'); plt.savefig('figures/Figure4_trials_and_synthesis.eps'); plt.close()
print('done')
