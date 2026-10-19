import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
import numpy as np, pandas as pd
from scipy import stats
from scipy.special import gammaln, logsumexp
T=pd.read_csv(f'{GW}/p3/trial_counts_from_sharifan_supp.tsv',sep='\t')
def lc(n,k): return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
g=np.linspace(np.log(0.05),np.log(20),4001); ll=np.zeros_like(g)
for _,t in T.iterrows():
    m=t.e1+t.e0; xs=np.arange(max(0,m-t.n0),min(m,t.n1)+1)
    ll+=lc(t.n1,t.e1)+lc(t.n0,t.e0)+t.e1*g-logsumexp((lc(t.n1,xs)+lc(t.n0,m-xs))[None,:]+np.outer(g,xs),axis=1)
ll-=ll.max()
def post(mu0,s0):
    lp=ll+stats.norm.logpdf(g,mu0,s0); p=np.exp(lp-lp.max()); p/=p.sum(); c=np.cumsum(p)
    return np.exp(g[np.searchsorted(c,0.5)]),np.exp(g[np.searchsorted(c,0.025)]),np.exp(g[np.searchsorted(c,0.975)]),p[g>0].sum(),p[g>np.log(1.2)].sum()
P=pd.read_csv(f'{GW}/p4/pooled_priors.tsv',sep='\t')
tot=P[P.prior=='total pooled'].iloc[0]; dir_=P[P.prior=='direct pooled'].iloc[0]
rows=[]
for k in [0.5,1,1.5,2,2.5,3]:
    for lab,(b,s) in [('face value, pooled PRACTICAL+FinnGen',(tot.beta,tot.se)),('PSA channel removed, pooled',(dir_.beta,dir_.se)),
                      ('face value, PRACTICAL only (sensitivity)',(0.1762,0.0511)),('PSA channel removed, PRACTICAL only (sensitivity)',(0.0397,0.1139))]:
        m,lo,hi,p1,p12=post(k*b,k*s); rows.append(dict(shift_SD=k,prior=lab,post_OR=m,lo=lo,hi=hi,P_gt1=p1,P_gt1_2=p12))
R=pd.DataFrame(rows); R.to_csv(f'{GW}/p4/synthesis_pooled.tsv',sep='\t',index=False)
pd.set_option('display.width',200); print(R[R.shift_SD.isin([1,2])].round(3).to_string(index=False))
z=(tot.beta-np.log(0.85))/np.sqrt(tot.se**2+((np.log(1.62)-np.log(0.44))/3.92)**2)
print('conflict, pooled face value vs RCT-only trials: z=%.2f p=%.3f'%(z,2*stats.norm.sf(abs(z))))
