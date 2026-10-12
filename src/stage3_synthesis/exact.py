import numpy as np,pandas as pd
from scipy import stats
from scipy.special import gammaln, logsumexp
T=[('Barnouin 2021',0,42,1,41),('Brock 2015',0,358,1,357),('Hackett 2013',0,92,1,98),('Marks 2006',2,21,4,19),('Wittert 2021',4,504,5,503),
   ('Rastrelli 2022',2,39,2,38),('Legros 2009',1,237,0,79),('Bhasin 2023 TRAVERSE',12,2601,11,2603),('Snyder 2016',2,395,1,395),('Basaria 2010',1,106,0,103)]
df=pd.DataFrame(T,columns=['trial','e1','n1','e0','n0']); df.to_csv('trial_counts_from_sharifan_supp.tsv',sep='\t',index=False)
def lc(n,k): return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
g=np.linspace(np.log(0.05),np.log(20),4001)
ll=np.zeros_like(g)
for _,r in df.iterrows():
    m=r.e1+r.e0
    xs=np.arange(max(0,m-r.n0),min(m,r.n1)+1)
    base=lc(r.n1,xs)+lc(r.n0,m-xs)
    num=lc(r.n1,r.e1)+lc(r.n0,r.e0)
    ll+=num+r.e1*g-logsumexp(base[None,:]+np.outer(g,xs),axis=1)
ll-=ll.max()
mle=g[ll.argmax()]; ci=g[ll>=-1.92]
print('exact conditional MLE OR %.3f, profile 95%% CI %.3f to %.3f, events %d/%d vs %d/%d'%(np.exp(mle),np.exp(ci.min()),np.exp(ci.max()),df.e1.sum(),df.n1.sum(),df.e0.sum(),df.n0.sum()))
dg=g[1]-g[0]
def post(mu0,s0):
    lp=ll+stats.norm.logpdf(g,mu0,s0); p=np.exp(lp-lp.max()); p/=p.sum()
    c=np.cumsum(p); med=g[np.searchsorted(c,0.5)]; lo=g[np.searchsorted(c,0.025)]; hi=g[np.searchsorted(c,0.975)]
    return np.exp(med),np.exp(lo),np.exp(hi),p[g>0].sum(),p[g>np.log(1.2)].sum()
rows=[]
for lab,(m0,s0) in [('vague',(0,10)),('face value MR, 1 SD',(0.1762,0.0511)),('face value MR, 2 SD',(0.3524,0.1022)),('PSA adjusted MR, 1 SD',(0.0397,0.1139)),('PSA adjusted MR, 2 SD',(0.0794,0.2278))]:
    r=post(m0,s0); rows.append(dict(prior=lab,post_OR=r[0],lo=r[1],hi=r[2],P_OR_gt1=r[3],P_OR_gt1_2=r[4]))
R=pd.DataFrame(rows); print(R.round(3).to_string(index=False)); R.to_csv('exact_posteriors.tsv',sep='\t',index=False)
# reverse Bayes with exact likelihood: prior SD equal to likelihood based SD
sd_lik=(ci.max()-ci.min())/3.92
for k in [0.5,1,2]:
    s0=sd_lik*k; lo,hi=-2,5
    for _ in range(60):
        mid=(lo+hi)/2
        if post(mid,s0)[3]>=0.95: hi=mid
        else: lo=mid
    print('reverse Bayes exact: prior SD %.3f (x%.1f of data) -> minimum prior median OR %.2f'%(s0,k,np.exp(hi)))
