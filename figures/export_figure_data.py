"""Writes the exact plotted values behind every figure panel to figures/data/,
one TSV per figure, so each image can be audited and redrawn without rerunning
the pipeline. Values are read from results/, never retyped.
"""
import pandas as pd, numpy as np, os
R='results/'; OUT='figures/data'; os.makedirs(OUT,exist_ok=True)
def mr(f,a,m):
    d=pd.read_csv(R+f,sep='\t'); r=d[(d.analysis==a)&(d.method==m)].iloc[0]; return r.beta,r.lo,r.hi,r.p,int(r.nsnp)

# Figure 1 is a conceptual schematic: its edges cite estimates shown elsewhere.
rows=[('bioavailable testosterone -> prostate AR signaling / PSA','raises PSA','F1, Figure 2'),
      ('PSA -> biopsy -> diagnosed prostate cancer','detection channel','F9, Figure 4A'),
      ('BMI -> PSA','lowers PSA','F9a, Figure 3'),
      ('androgen signaling -> tumour initiation -> lethal cancer','not demonstrated','L2, untested')]
pd.DataFrame(rows,columns=['edge','label','evidence']).to_csv(f'{OUT}/Figure1_conceptual_model_edges.tsv',sep='\t',index=False)

# Figure 2: saturation and free hormone test
MX=pd.read_csv(R+'stage7_matrix/estimator_matrix.tsv',sep='\t')
f2=[]
for lab,a,e in [('Bioavailable T, IVW','biot->psa','IVW'),
                ('Bioavailable T, weighted median','biot->psa','Weighted median'),
                ('Bioavailable T, weighted mode','biot->psa','Weighted mode'),
                ('Bioavailable T, IVW Q outliers removed','biot->psa','IVW, Q outliers removed'),
                ('Bioavailable T, MR cML (DP SE)','biot->psa','MR cML BIC (DP SE)'),
                ('Bioavailable T, MR RAPS (bootstrap SE)','biot->psa','MR RAPS Huber'),
                ('Bioavailable T, MR Egger','biot->psa','MR Egger'),
                ('Total T, IVW','tt->psa','IVW'),
                ('SHBG, IVW','shbg->psa','IVW')]:
    r=MX[(MX.analysis==a)&(MX.estimator==e)].iloc[0]
    f2.append(dict(panel_row=lab,beta=r.beta,lo=r.lo,hi=r.hi,p=r.p,nsnp=int(r.nsnp),
                   source='stage7_matrix/estimator_matrix.tsv'))
pd.DataFrame(f2).to_csv(f'{OUT}/Figure2_saturation_test_data.tsv',sep='\t',index=False)

# Figure 3: detection signature, three exposures x three outcomes
U=pd.read_csv(R+'stage2/uni_metabolic.tsv',sep='\t'); UP=pd.read_csv(R+'stage2/uni_metabolic_psa.tsv',sep='\t')
f3=[]
for lab,k in [('Bioavailable testosterone','biot'),('BMI','bmi'),('IGF-1','igf')]:
    for panel,outc in [('Effect on PSA (log units)','psa'),('Prostate cancer, PRACTICAL (OR)','prac'),('Prostate cancer, FinnGen (OR)','fg')]:
        src='stage1/mr_results.tsv' if k=='biot' else ('stage2/uni_metabolic_psa.tsv' if outc=='psa' else 'stage2/uni_metabolic.tsv')
        b,lo,hi,p,n=mr(src,f'{k}->{outc}','IVW')
        if outc!='psa': b,lo,hi=np.exp([b,lo,hi])
        f3.append(dict(exposure=lab,panel=panel,estimate=b,lo=lo,hi=hi,p=p,nsnp=n,scale='log units' if outc=='psa' else 'odds ratio',source=src))
pd.DataFrame(f3).to_csv(f'{OUT}/Figure3_detection_signature_data.tsv',sep='\t',index=False)

# Figure 4A: trial inputs (published, not re-derived here)
pd.DataFrame([dict(outcome='Prostate biopsies',RR=2.38,lo=1.01,hi=5.60,source='Zeng 2022, randomised trials, doi:10.3389/fendo.2022.1009900'),
              dict(outcome='Prostate cancer',RR=0.85,lo=0.44,hi=1.62,source='Zeng 2022, randomised trials only sensitivity analysis')]
            ).to_csv(f'{OUT}/Figure4A_trial_inputs.tsv',sep='\t',index=False)
# Figure 4B: posterior curves as plotted
S=pd.read_csv(R+'stage4/synthesis_pooled.tsv',sep='\t')
S[S.prior.isin(['face value, pooled PRACTICAL+FinnGen','PSA channel removed, pooled'])].to_csv(f'{OUT}/Figure4B_posterior_curves.tsv',sep='\t',index=False)
for f in sorted(os.listdir(OUT)): print('wrote', f'{OUT}/{f}')
