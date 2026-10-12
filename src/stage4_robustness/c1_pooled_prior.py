import pandas as pd, numpy as np, sys
from scipy import stats
R=pd.read_csv('/tmp/gw/out/mr_results.tsv',sep='\t')
def g(a,m): r=R[(R.analysis==a)&(R.method==m)].iloc[0]; return r.beta,r.se
bp,sp=g('biot->prac','IVW'); bf,sf=g('biot->fg','IVW')
w=[1/sp**2,1/sf**2]; b=(bp*w[0]+bf*w[1])/sum(w); se=np.sqrt(1/sum(w))
Q=w[0]*(bp-b)**2+w[1]*(bf-b)**2; pq=stats.chi2.sf(Q,1); I2=max(0,(Q-1)/Q)*100 if Q>0 else 0
print('PRACTICAL  %.4f (%.4f)  OR %.3f'%(bp,sp,np.exp(bp)))
print('FinnGen    %.4f (%.4f)  OR %.3f'%(bf,sf,np.exp(bf)))
print('POOLED     %.4f (%.4f)  OR %.3f (%.3f-%.3f) p=%.4f  Q=%.2f p=%.3f I2=%.0f%%'%(b,se,np.exp(b),np.exp(b-1.96*se),np.exp(b+1.96*se),2*stats.norm.sf(abs(b/se)),Q,pq,I2))
# detection-adjusted pooled: apply the same proportional attenuation observed in MVMR
MV=pd.read_csv('/tmp/gw/out/mvmr_results.tsv',sep='\t')
dp=MV[(MV.model=='biot+psa->prac')&(MV.exposure=='biot')].iloc[0]
df_=MV[(MV.model=='biot+psa->fg')&(MV.exposure=='biot')].iloc[0]
w2=[1/dp.se**2,1/df_.se**2]; bd=(dp.beta*w2[0]+df_.beta*w2[1])/sum(w2); sd_=np.sqrt(1/sum(w2))
Qd=w2[0]*(dp.beta-bd)**2+w2[1]*(df_.beta-bd)**2
print('DIRECT pooled %.4f (%.4f) OR %.3f (%.3f-%.3f) Q p=%.3f'%(bd,sd_,np.exp(bd),np.exp(bd-1.96*sd_),np.exp(bd+1.96*sd_),stats.chi2.sf(Qd,1)))
pd.DataFrame([dict(prior='total pooled',beta=b,se=se,Q=Q,Q_p=pq),dict(prior='direct pooled',beta=bd,se=sd_,Q=Qd,Q_p=stats.chi2.sf(Qd,1))]).to_csv('/tmp/gw/p4/pooled_priors.tsv',sep='\t',index=False)
# I1 heterogeneity for metabolic exposures
U=pd.read_csv('/tmp/gw/p2/uni_metabolic.tsv',sep='\t')
rows=[]
for t,lab in [('bmi','BMI'),('igf','IGF-1'),('fi','Fasting insulin')]:
    a=U[(U.analysis==f'{t}->prac')&(U.method=='IVW')].iloc[0]; c=U[(U.analysis==f'{t}->fg')&(U.method=='IVW')].iloc[0]
    z=(a.beta-c.beta)/np.sqrt(a.se**2+c.se**2); rows.append(dict(exposure=lab,prac_OR=np.exp(a.beta),fg_OR=np.exp(c.beta),het_z=z,het_p=2*stats.norm.sf(abs(z))))
rows.append(dict(exposure='Bioavailable testosterone',prac_OR=np.exp(bp),fg_OR=np.exp(bf),het_z=(bp-bf)/np.sqrt(sp**2+sf**2),het_p=pq))
H=pd.DataFrame(rows); H.to_csv('/tmp/gw/p4/cohort_heterogeneity.tsv',sep='\t',index=False); print(H.round(3).to_string(index=False))
