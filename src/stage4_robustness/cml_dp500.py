import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
import numpy as np,pandas as pd,sys
from scipy import stats
sys.path.insert(0,GW); import mr
exec(open(f'{GW}/robust2.py').read().split('# ---------- MR RAPS')[0].split("rng=np.random.default_rng(20261011)")[1])
rng=np.random.default_rng(7)
def cml_fast(bx,by,sx,sy,nmin,Kmax):
    best=None; init=np.sum(bx*by/sy**2)/np.sum(bx**2/sy**2)
    for K in range(0,Kmax+1):
        f=min((cml_fit(bx,by,sx,sy,K,s,iters=100) for s in [init,init-0.2,init+0.2]),key=lambda z:z[2])
        bic=2*f[2]+np.log(nmin)*K
        if best is None or bic<best[1]: best=(f[0],bic)
    return best[0]
D=mr.load(); out=[]
for e,o,nmin,tag in [('biot','psa',178782,'psa'),('biot','prac',140254,'prac'),('biot','fg',178782,'fg')]:
    inst=set(pd.read_csv(f'{GW}/{e}_clump.clumped',sep=r'\s+').SNP)
    m=mr.harmonise(D[e][D[e].SNP.isin(inst)],D[o]); bx,by,sx,sy=m.b_x.values,m.b_y.values,m.se_x.values,m.se_y.values
    est=cml_fast(bx,by,sx,sy,nmin,25)
    bs=[cml_fast(bx+rng.normal(0,sx),by+rng.normal(0,sy),sx,sy,nmin,25) for _ in range(500)]
    se=np.std(bs); out.append(dict(analysis=f'{e}->{o}',method='MR cML BIC, data perturbation SE (B=500)',beta=est,se=se,lo=est-1.96*se,hi=est+1.96*se,p=2*stats.norm.sf(abs(est/se))))
    print(out[-1],flush=True)
pd.DataFrame(out).to_csv(f'{GW}/p4/cml_dp_500.tsv',sep='\t',index=False)
