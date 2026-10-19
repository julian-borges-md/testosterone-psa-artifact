import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
import pandas as pd,numpy as np,mr,subprocess
from scipy import stats
from scipy.optimize import minimize_scalar
D=mr.load()
def joint(names,out):
    # joint instrument set: union of genome-wide hits, re-clumped on min p
    frames=[]
    for t in names:
        c=pd.read_csv(f'{t}_clump.clumped',sep=r'\s+')[['SNP','P']]; frames.append(c)
    u=pd.concat(frames).groupby('SNP',as_index=False).P.min()
    u.to_csv('j_in.txt',sep='\t',index=False)
    subprocess.run([f'{TOOLS}/plink','--bfile',f'{TOOLS}/EUR','--clump','j_in.txt','--clump-p1','5e-8','--clump-p2','1','--clump-r2','0.001','--clump-kb','10000','--out','j','--silent','--memory','2500'],check=True)
    snps=set(pd.read_csv('j.clumped',sep=r'\s+').SNP)
    base=D[names[0]][D[names[0]].SNP.isin(snps)]
    m=mr.harmonise(base,D[out]).rename(columns={'b_x':'bx0','se_x':'sx0','b_y':'by','se_y':'sy'})[['SNP','ea_x','oa_x','bx0','sx0','by','sy']]
    for i,t in enumerate(names[1:],1):
        h=mr.harmonise(base,D[t]).rename(columns={'b_y':f'bx{i}','se_y':f'sx{i}'})[['SNP',f'bx{i}',f'sx{i}']]
        m=m.merge(h,on='SNP')
    return m,len(snps)
def mvmr(m,names):
    k=len(m); X=np.column_stack([m[f'bx{i}'] for i in range(len(names))]); w=1/m.sy.values**2
    XtW=X.T*w; b=np.linalg.solve(XtW@X,XtW@m.by.values); res=m.by.values-X@b
    phi=max(1,np.sum(w*res**2)/(k-len(names))); se=np.sqrt(np.diag(np.linalg.inv(XtW@X))*phi)
    # conditional F (two exposures)
    cf=[]
    for j in range(2):
        o=1-j; bj,sj,bo,so=m[f'bx{j}'].values,m[f'sx{j}'].values,m[f'bx{o}'].values,m[f'sx{o}'].values
        f=lambda d: np.sum((bj-d*bo)**2/(sj**2+d**2*so**2))
        d=minimize_scalar(f).x; cf.append(f(d)/(k-1))
    return b,se,cf,k
rows=[]
for names,out in [(['biot','psa'],'prac'),(['biot','psa'],'fg'),(['biot','shbg'],'psa'),(['biot','shbg'],'prac')]:
    m,nj=joint(names,out); b,se,cf,k=mvmr(m,names)
    for i,n in enumerate(names):
        lo,hi=b[i]-1.96*se[i],b[i]+1.96*se[i]; p=2*stats.norm.sf(abs(b[i]/se[i]))
        rows.append(dict(model='+'.join(names)+'->'+out,exposure=n,nsnp=k,beta=b[i],se=se[i],lo=lo,hi=hi,p=p,condF=cf[i]))
r=pd.DataFrame(rows); r.to_csv('out/mvmr_results.tsv',sep='\t',index=False)
pd.set_option('display.width',200); print(r.round(4).to_string(index=False))
