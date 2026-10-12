import numpy as np, pandas as pd, sys
from scipy import stats, optimize
sys.path.insert(0,'/tmp/gw'); import mr
rng=np.random.default_rng(20261011)

# ---------- MR cML (Xue, Shen, Pan 2021), BIC selection and model averaging ----------
def cml_fit(bx,by,sx,sy,K,theta0,iters=200):
    th=theta0; n=len(bx); r=np.zeros(n)
    for _ in range(iters):
        bj=(bx/sx**2+th*(by-r)/sy**2)/(1/sx**2+th**2/sy**2)
        res=(by-th*bj)**2/sy**2
        inv=np.argsort(-res)[:K] if K>0 else np.array([],int)
        r=np.zeros(n); r[inv]=by[inv]-th*bj[inv]
        v=np.ones(n,bool); v[inv]=False
        th_new=np.sum(bj[v]*by[v]/sy[v]**2)/np.sum(bj[v]**2/sy[v]**2)
        if abs(th_new-th)<1e-10: th=th_new; break
        th=th_new
    bj=(bx/sx**2+th*(by-r)/sy**2)/(1/sx**2+th**2/sy**2)
    nll=0.5*np.sum((bx-bj)**2/sx**2)+0.5*np.sum((by-th*bj-r)**2/sy**2)
    v=r==0
    se=1/np.sqrt(np.sum(bx[v]**2/(sy[v]**2+th**2*sx[v]**2)))
    return th,se,nll,v
def cml(bx,by,sx,sy,nmin,starts=10):
    n=len(bx); fits=[]
    init=np.sum(bx*by/sy**2)/np.sum(bx**2/sy**2)
    for K in range(0,n-2):
        best=None
        for s in np.concatenate([[init],rng.uniform(init-0.5,init+0.5,starts)]):
            f=cml_fit(bx,by,sx,sy,K,s)
            if best is None or f[2]<best[2]: best=f
        bic=2*best[2]+np.log(nmin)*K
        fits.append((K,best[0],best[1],bic,best[3]))
        if K>n//2: break
    F=pd.DataFrame(fits,columns=['K','theta','se','bic','valid'])
    b=F.loc[F.bic.idxmin()]
    w=np.exp(-(F.bic-F.bic.min())/2); w/=w.sum()
    ma=np.sum(w*F.theta); mase=np.sqrt(np.sum(w*np.sqrt(F.se**2+(F.theta-ma)**2))**2)
    return dict(K_bic=int(b.K),theta_bic=b.theta,se_bic=b.se,theta_ma=ma,se_ma=mase)

# ---------- MR RAPS with Huber loss (Zhao et al. 2020), bootstrap SE ----------
def raps_est(bx,by,sx,sy,c=1.345):
    def obj(p):
        th,lt=p; tau2=np.exp(lt)
        t=(by-th*bx)/np.sqrt(sy**2+th**2*sx**2+tau2)
        a=np.abs(t); rho=np.where(a<=c,0.5*t**2,c*a-0.5*c**2)
        return np.sum(rho)+0.5*np.sum(np.log(sy**2+th**2*sx**2+tau2))*0.0+0.0
    def obj_full(p):
        th,lt=p; tau2=np.exp(lt); V=sy**2+th**2*sx**2+tau2
        t=(by-th*bx)/np.sqrt(V); a=np.abs(t); rho=np.where(a<=c,0.5*t**2,c*a-0.5*c**2)
        return np.sum(rho)+0.5*np.sum(np.log(V))
    init=np.sum(bx*by/sy**2)/np.sum(bx**2/sy**2)
    r=optimize.minimize(obj_full,[init,np.log(np.median(sy**2))],method='Nelder-Mead',options={'xatol':1e-8,'fatol':1e-10,'maxiter':4000})
    return r.x[0],np.exp(r.x[1])
def raps(bx,by,sx,sy,B=300):
    th,tau2=raps_est(bx,by,sx,sy)
    bs=[]
    for _ in range(B):
        i=rng.integers(0,len(bx),len(bx))
        bs.append(raps_est(bx[i],by[i],sx[i],sy[i])[0])
    return th,np.std(bs),tau2

D=mr.load(); inst={t:set(pd.read_csv(f'/tmp/gw/{t}_clump.clumped',sep=r'\s+').SNP) for t in ['biot','shbg','tt','psa']}
rows=[]
for e,o,nmin in [('biot','psa',178782),('biot','prac',140254),('biot','fg',178782),('shbg','psa',180726),('tt','psa',194453),('psa','prac',140254)]:
    ex=D[e][D[e].SNP.isin(inst[e])]; m=mr.harmonise(ex,D[o])
    bx,by,sx,sy=m.b_x.values,m.b_y.values,m.se_x.values,m.se_y.values
    ivw=mr.ivw(bx,by,sy)
    c=cml(bx,by,sx,sy,nmin)
    rt,rse,tau2=raps(bx,by,sx,sy)
    for meth,b,s in [('IVW (reference)',ivw[0],ivw[1]),('MR cML BIC',c['theta_bic'],c['se_bic']),('MR cML model averaged',c['theta_ma'],c['se_ma']),('MR RAPS Huber',rt,rse)]:
        rows.append(dict(analysis=f'{e}->{o}',method=meth,nsnp=len(m),beta=b,se=s,lo=b-1.96*s,hi=b+1.96*s,p=2*stats.norm.sf(abs(b/s)),K_invalid=c['K_bic'] if 'cML' in meth else np.nan,tau2=tau2 if 'RAPS' in meth else np.nan))
    print(e,o,'done',flush=True)
R=pd.DataFrame(rows); R.to_csv('/tmp/gw/p2/robust_cml_raps.tsv',sep='\t',index=False)
pd.set_option('display.width',220); print(R.round(4).to_string(index=False))
