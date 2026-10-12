import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(20261009)
bim=pd.read_csv('/tmp/gw/union_snps.tsv',sep='\t',dtype={'chr':str})
def std(df,ea,oa,eaf,b,se,p):
    d=pd.DataFrame({'SNP':df.SNP,'ea':df[ea].astype(str).str.upper(),'oa':df[oa].astype(str).str.upper(),
       'eaf':pd.to_numeric(df[eaf],errors='coerce'),'b':pd.to_numeric(df[b],errors='coerce'),'se':pd.to_numeric(df[se],errors='coerce'),'p':pd.to_numeric(df[p],errors='coerce')})
    d=d[d.ea.str.len().eq(1)&d.oa.str.len().eq(1)].dropna(subset=['b','se'])
    return d.drop_duplicates('SNP')
def load():
    D={}
    for t,n in [('biot',178782),('shbg',180726),('tt',194453)]:
        x=pd.read_csv(f'{t}_u.tsv',sep='\t'); D[t]=std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value'); D[t]['n']=n
    x=pd.read_csv('psa_u.tsv',sep='\t'); d=std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')
    # PSA rows matched on position: keep rows whose alleles match reference
    ref=bim.set_index('rsid')[['a1','a2']]
    d=d.join(ref,on='SNP'); d=d[((d.ea==d.a1)&(d.oa==d.a2))|((d.ea==d.a2)&(d.oa==d.a1))].drop(columns=['a1','a2'])
    d['n']=pd.to_numeric(x.drop_duplicates('SNP').set_index('SNP').loc[d.SNP,'n'].values,errors='coerce'); D['psa']=d
    x=pd.read_csv('prac_u.tsv',sep='\t'); d=std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')
    d=d.join(ref,on='SNP'); d=d[((d.ea==d.a1)&(d.oa==d.a2))|((d.ea==d.a2)&(d.oa==d.a1))].drop(columns=['a1','a2']); d['n']=140254; D['prac']=d
    x=pd.read_csv('fg_u.tsv',sep='\t'); D['fg']=std(x,'alt','ref','af_alt','beta','sebeta','pval'); D['fg']['n']=np.nan
    return D
PAL={('A','T'),('T','A'),('C','G'),('G','C')}
def harmonise(exp,out,log=None):
    m=exp.merge(out,on='SNP',suffixes=('_x','_y')); n0=len(m)
    same=(m.ea_x==m.ea_y)&(m.oa_x==m.oa_y); flip=(m.ea_x==m.oa_y)&(m.oa_x==m.ea_y)
    m=m[same|flip].copy(); fl=flip[same|flip]
    m.loc[fl,'b_y']*=-1; m.loc[fl,'eaf_y']=1-m.loc[fl,'eaf_y']
    n1=len(m)
    pal=[(a,b) in PAL for a,b in zip(m.ea_x,m.oa_x)]; maf=np.minimum(m.eaf_x,1-m.eaf_x)
    m=m[~(np.array(pal)&(maf>0.42))]; n2=len(m)
    disc=(m.eaf_y.notna())&((m.eaf_x-m.eaf_y).abs()>0.2); m=m[~disc]; n3=len(m)
    if log is not None: log.update(dict(merged=n0,allele_ok=n1,pal_ok=n2,eaf_ok=n3))
    return m
def ivw(bx,by,se):
    w=1/se**2; r=by/bx; sw=se/np.abs(bx)
    W=1/sw**2; b=np.sum(W*r)/np.sum(W)
    Q=np.sum(W*(r-b)**2); k=len(r); phi=max(1,Q/(k-1))
    s=np.sqrt(1/np.sum(W))*np.sqrt(phi); return b,s,Q,k
def egger(bx,by,se):
    s=np.sign(bx); bx=bx*s; by=by*s; w=1/se**2
    X=np.column_stack([np.ones_like(bx),bx]); W=np.diag(w)
    XtWX=X.T@W@X; beta=np.linalg.solve(XtWX,X.T@W@by)
    res=by-X@beta; k=len(bx); sig2=max(1,np.sum(w*res**2)/(k-2))
    cov=np.linalg.inv(XtWX)*sig2; return beta[1],np.sqrt(cov[1,1]),beta[0],np.sqrt(cov[0,0])
def wmed_est(r,w):
    o=np.argsort(r); r=r[o]; w=w[o]/w.sum(); c=np.cumsum(w)-0.5*w
    i=np.searchsorted(c,0.5); 
    if i==0: return r[0]
    return r[i-1]+(r[i]-r[i-1])*(0.5-c[i-1])/(c[i]-c[i-1])
def wmed(bx,by,sey,sex,B=1000):
    r=by/bx; w=1/(sey/np.abs(bx))**2; e=wmed_est(r,w)
    bs=[wmed_est((by+rng.normal(0,sey))/(bx+rng.normal(0,sex)),w) for _ in range(B)]
    return e,np.std(bs)
def wmode_est(r,se_r):
    w=(se_r**-2)/np.sum(se_r**-2); s=0.9*min(np.std(r),stats.median_abs_deviation(r,scale='normal'))/len(r)**0.2
    g=np.linspace(r.min(),r.max(),1000); dens=(w[None,:]*np.exp(-0.5*((g[:,None]-r[None,:])/s)**2)).sum(1); return g[dens.argmax()]
def wmode(bx,by,sey,sex,B=200):
    r=by/bx; se_r=sey/np.abs(bx); e=wmode_est(r,se_r)
    bs=[wmode_est((by+rng.normal(0,sey))/(bx+rng.normal(0,sex)),se_r) for _ in range(B)]
    return e,np.std(bs)
def ci(b,s): return b-1.96*s,b+1.96*s,2*stats.norm.sf(abs(b/s))
def run(m,label):
    bx,by,sx,sy=m.b_x.values,m.b_y.values,m.se_x.values,m.se_y.values
    F=np.mean((bx/sx)**2)
    b,s,Q,k=ivw(bx,by,sy); Qp=stats.chi2.sf(Q,k-1)
    rows=[('IVW',b,s)]
    eb,es,ei,eis=egger(bx,by,sy); rows.append(('MR Egger',eb,es))
    wb,ws=wmed(bx,by,sy,sx); rows.append(('Weighted median',wb,ws))
    mb,ms=wmode(bx,by,sy,sx); rows.append(('Weighted mode',mb,ms))
    # Q outlier removal (Bonferroni), then IVW
    W=(bx/sy)**2; q=W*(by/bx-b)**2; keep=stats.chi2.sf(q,1)>0.05/k
    if keep.sum()<k and keep.sum()>2:
        ob,os_,_,ok=ivw(bx[keep],by[keep],sy[keep]); rows.append((f'IVW outliers removed ({k-keep.sum()})',ob,os_))
    out=[]
    for meth,bb,ss in rows:
        lo,hi,p=ci(bb,ss); out.append(dict(analysis=label,method=meth,nsnp=k,beta=bb,se=ss,lo=lo,hi=hi,p=p))
    egp=2*stats.norm.sf(abs(ei/eis))
    meta=dict(analysis=label,nsnp=k,meanF=F,Q=Q,Qp=Qp,egger_intercept=ei,egger_int_p=egp)
    return out,meta
