import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import NullFormatter, NullLocator
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
C={'biot':'#1f4e79','tt':'#7f7f7f','shbg':'#c55a11','bmi':'#548235','igf':'#7030a0'}
R=pd.read_csv('results/stage1/mr_results.tsv',sep='\t'); RB=pd.read_csv('results/stage2/robust_cml_raps.tsv',sep='\t'); DP=pd.read_csv('results/stage2/cml_dp.tsv',sep='\t')
U=pd.read_csv('results/stage2/uni_metabolic.tsv',sep='\t'); UP=pd.read_csv('results/stage2/uni_metabolic_psa.tsv',sep='\t')
S=pd.read_csv('results/stage4/synthesis_pooled.tsv',sep='\t')
MX=pd.read_csv('results/stage7_matrix/estimator_matrix.tsv',sep='\t')
def mx(a,e):
    r=MX[(MX.analysis==a)&(MX.estimator==e)].iloc[0]; return r.beta,r.lo,r.hi
def get(df,a,m): r=df[(df.analysis==a)&(df.method==m)].iloc[0]; return r.beta,r.lo,r.hi

# ---------- Figure 1: conceptual model ----------
fig,ax=plt.subplots(figsize=(7.6,3.8)); ax.axis('off'); ax.set_xlim(0,10.4); ax.set_ylim(0,5)
def box(x,y,t,fc): ax.add_patch(FancyBboxPatch((x-0.95,y-0.38),1.9,0.76,boxstyle='round,pad=0.05',fc=fc,ec='black',lw=0.8)); ax.text(x,y,t,ha='center',va='center',fontsize=8)
def arr(a,b,t='',c='black',ls='-',dy=0.15):
    ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',color=c,lw=1.2,ls=ls))
    if t: ax.text((a[0]+b[0])/2,(a[1]+b[1])/2+dy,t,ha='center',fontsize=7,color=c)
box(1.1,3.5,'Bioavailable\ntestosterone','#dbe5f1'); box(3.6,3.5,'Prostate AR\nsignaling','#dbe5f1')
box(6.4,4.3,'PSA','#fde9d9'); box(9.2,4.3,'Biopsy','#fde9d9'); box(9.2,2.6,'Diagnosed\nprostate cancer','#f2f2f2')
box(6.4,1.4,'Tumor initiation\nor progression','#e2efda'); box(9.2,1.0,'Lethal prostate\ncancer','#e2efda')
box(1.1,1.6,'BMI\n(hemodilution)','#ebf1de')
arr((2.05,3.5),(2.65,3.5)); arr((4.55,3.75),(5.45,4.2)); ax.text(4.75,4.35,'raises PSA\n(saturation falsified)',ha='center',fontsize=7,color='#1f4e79')
arr((7.35,4.3),(8.25,4.3),c='#c00000'); arr((9.2,3.92),(9.2,2.98),c='#c00000'); ax.text(9.35,3.45,'detection',fontsize=7,color='#c00000')
arr((4.4,3.12),(5.6,1.75),c='#548235',ls='--'); ax.text(4.55,2.3,'biology?',fontsize=7,color='#548235')
arr((7.35,1.6),(8.25,2.35),c='#548235',ls='--'); arr((7.35,1.2),(8.25,1.0),c='#548235',ls='--')
arr((2.05,1.75),(5.45,4.0),c='#7f7f7f'); ax.text(2.9,2.75,'lowers PSA',fontsize=7,color='#7f7f7f',rotation=30)
ax.text(5.2,0.1,'Red path: PSA driven detection (supported). Dashed green path: tumor biology (not demonstrated for lethal disease).',ha='center',fontsize=7,style='italic')
plt.savefig('figures/Figure1_conceptual_model.png',bbox_inches='tight'); plt.savefig('figures/Figure1_conceptual_model.pdf',bbox_inches='tight'); plt.close()

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
fig,ax=plt.subplots(figsize=(6.6,4.0)); y=np.arange(len(rows))[::-1]
for yi,(lab,b,lo,hi,c) in zip(y,rows):
    ax.errorbar(b,yi,xerr=[[b-lo],[hi-b]],fmt='o',color=c,ms=5,capsize=2,lw=1.2)
ax.axvline(0,color='black',lw=0.8,ls=':')
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.set_xlabel('Effect on PSA (natural log units per SD of exposure, 95% CI)')
ax.text(0.003,len(rows)-0.4,'Saturation predicts zero',fontsize=7,style='italic')
plt.tight_layout(); plt.savefig('figures/Figure2_saturation_test.png'); plt.savefig('figures/Figure2_saturation_test.pdf'); plt.close()

# ---------- Figure 3: detection signature ----------
ex=[('Bioavailable\ntestosterone','biot',C['biot']),('BMI','bmi',C['bmi']),('IGF-1','igf',C['igf'])]
fig,axs=plt.subplots(1,3,figsize=(7.2,2.6),sharey=True)
titles=['Effect on PSA\n(log units)','Prostate cancer, PRACTICAL\n(OR; heavily screened)','Prostate cancer, FinnGen\n(OR; registry)']
for j,ax in enumerate(axs):
    for i,(lab,k,c) in enumerate(ex):
        if j==0: src=R if k=='biot' else UP; a=f'{k}->psa'
        else:
            o='prac' if j==1 else 'fg'; src=R if k=='biot' else U; a=f'{k}->{o}'
        b,lo,hi=get(src,a,'IVW')
        if j>0: b,lo,hi=np.exp([b,lo,hi])
        yi=2-i; ax.errorbar(b,yi,xerr=[[b-lo],[hi-b]],fmt='o',color=c,ms=5,capsize=2,lw=1.2)
    ax.axvline(0 if j==0 else 1,color='black',lw=0.8,ls=':'); ax.set_title(titles[j],fontsize=8)
axs[0].set_yticks([2,1,0]); axs[0].set_yticklabels([e[0] for e in ex])
plt.tight_layout(); plt.savefig('figures/Figure3_detection_signature.png'); plt.savefig('figures/Figure3_detection_signature.pdf'); plt.close()

# ---------- Figure 4: trials and synthesis ----------
fig,axs=plt.subplots(1,2,figsize=(7.2,2.9),gridspec_kw={'width_ratios':[1,1.3]})
ax=axs[0]; tr=[('Prostate biopsies',2.38,1.01,5.60,'#c00000'),('Prostate cancer',0.85,0.44,1.62,'#1f4e79')]
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
ax.set_yticks([1,0]); ax.set_yticklabels([t[0] for t in tr]); ax.set_ylim(-0.6,1.6); ax.set_xlabel('Relative risk, randomized trials (95% CI)'); ax.set_title('A. Testosterone therapy in trials',fontsize=8,loc='left')
ax=axs[1]; d=S
for pr,c,lab in [('face value, pooled PRACTICAL+FinnGen','#7f7f7f','Genetics at face value'),('PSA channel removed, pooled','#1f4e79','Genetics, PSA channel removed')]:
    q=d[d.prior==pr]; ax.plot(q.shift_SD,q.post_OR,'-o',color=c,ms=3,label=lab); ax.fill_between(q.shift_SD,q.lo,q.hi,color=c,alpha=0.15)
ax.axhline(1,color='black',lw=0.8,ls=':'); ax.axhline(1.2,color='#c00000',lw=0.7,ls='--'); ax.text(0.55,1.22,'meaningful harm (OR 1.2)',fontsize=7,color='#c00000',
        bbox=dict(facecolor='white',edgecolor='none',pad=1.0))
ax.set_xlabel('Assumed therapy shift (SD)'); ax.set_ylabel('Posterior OR, any prostate cancer'); ax.legend(fontsize=7,frameon=False,loc='upper left')
ax.set_title('B. Bayesian synthesis with 41 trials',fontsize=8,loc='left')
plt.tight_layout(); plt.savefig('figures/Figure4_trials_and_synthesis.png'); plt.savefig('figures/Figure4_trials_and_synthesis.pdf'); plt.close()
print('done')
