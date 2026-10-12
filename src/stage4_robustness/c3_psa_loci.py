import pandas as pd, numpy as np, sys, subprocess
sys.path.insert(0,'/tmp/gw'); import mr
from scipy import stats
from scipy.optimize import minimize_scalar
# GRCh37 1Mb windows around canonical prostate-cancer / PSA-biology loci
LOCI=[('KLK3','19',50300000,52400000),('MSMB','10',50500000,52600000),('TERT','5',250000,2300000),
      ('HNF1B','17',35000000,37100000),('SHBG','17',6500000,8600000)]
bim=pd.read_csv('/tmp/tools/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str}).set_index('rsid')
def flag(snps):
    r=bim.loc[snps]; m=pd.Series(False,index=snps)
    for g,c,a,b in LOCI: m |= ((r.chr==c)&(r.pos>=a)&(r.pos<=b)).values
    return m.values
D=mr.load()
clump={'biot':'/tmp/gw/biot_clump.clumped','psa':'/tmp/gw/psa_clump.clumped'}
def joint(names):
    u=pd.concat([pd.read_csv(clump[t],sep=r'\s+')[['SNP','P']] for t in names]).groupby('SNP',as_index=False).P.min()
    u.to_csv('/tmp/gw/p4/jin.txt',sep='\t',index=False)
    subprocess.run(['/tmp/tools/plink','--bfile','/tmp/tools/EUR','--clump','/tmp/gw/p4/jin.txt','--clump-p1','5e-8','--clump-p2','1','--clump-r2','0.001','--clump-kb','10000','--out','/tmp/gw/p4/jj','--silent','--memory','2000'],check=True)
    return set(pd.read_csv('/tmp/gw/p4/jj.clumped',sep=r'\s+').SNP)
def mvmr(m):
    X=np.column_stack([m.bx0,m.bx1]); w=1/m.sy.values**2; XtW=X.T*w
    b=np.linalg.solve(XtW@X,XtW@m.by.values); res=m.by.values-X@b; k=len(m)
    phi=max(1,np.sum(w*res**2)/(k-2)); se=np.sqrt(np.diag(np.linalg.inv(XtW@X))*phi)
    cf=[]
    for j,o in [(0,1),(1,0)]:
        bj,sj,bo,so=m[f'bx{j}'].values,m[f'sx{j}'].values,m[f'bx{o}'].values,m[f'sx{o}'].values
        f=lambda d: np.sum((bj-d*bo)**2/(sj**2+d**2*so**2)); d=minimize_scalar(f).x; cf.append(f(d)/(k-1))
    return b,se,cf,k
out=[]
for outc in ['prac','fg']:
    snps=joint(['biot','psa']); base=D['biot'][D['biot'].SNP.isin(snps)]
    m=mr.harmonise(base,D[outc]).rename(columns={'b_x':'bx0','se_x':'sx0','b_y':'by','se_y':'sy'})[['SNP','bx0','sx0','by','sy']]
    h=mr.harmonise(base,D['psa']).rename(columns={'b_y':'bx1','se_y':'sx1'})[['SNP','bx1','sx1']]
    m=m.merge(h,on='SNP')
    for lab,sub in [('all instruments',m),('excluding PSA biology and PCa loci',m[~flag(m.SNP.values)])]:
        b,se,cf,k=mvmr(sub)
        for i,e in enumerate(['bioavailable testosterone (direct)','PSA']):
            out.append(dict(outcome=outc,set=lab,exposure=e,nsnp=k,OR=np.exp(b[i]),lo=np.exp(b[i]-1.96*se[i]),hi=np.exp(b[i]+1.96*se[i]),p=2*stats.norm.sf(abs(b[i]/se[i])),condF=cf[i]))
        if outc=='prac':
            a=mr.ivw(*[mr.harmonise(D['biot'][D['biot'].SNP.isin(set(pd.read_csv(clump['biot'],sep=r'\s+').SNP))],D['psa'])[c].values for c in ['b_x','b_y','se_y']])[0]
            tot=mr.ivw(*[mr.harmonise(D['biot'][D['biot'].SNP.isin(set(pd.read_csv(clump['biot'],sep=r'\s+').SNP))],D['prac'])[c].values for c in ['b_x','b_y','se_y']])[0]
            print(f'  {lab}: PSA channel share = {a*b[1]/tot:.3f}')
O=pd.DataFrame(out); O.to_csv('/tmp/gw/p4/c3_psa_loci_sensitivity.tsv',sep='\t',index=False)
pd.set_option('display.width',200); print(O.round(3).to_string(index=False))
