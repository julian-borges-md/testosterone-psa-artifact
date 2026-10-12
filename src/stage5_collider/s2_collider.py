"""S2: is conditioning on PSA mediation or collider bias?
Negative-control logic: if conditioning on PSA were chiefly collider bias, it would
distort exposures regardless of whether they act through PSA. We test an exposure
with no PSA effect (IGF-1) and one with an opposite PSA effect (BMI).
Prespecified: if MVMR leaves IGF-1 unchanged but moves PSA-acting exposures, the
adjustment behaves as mediation, not as a generic collider distortion.
"""
import pandas as pd, numpy as np, sys, subprocess
sys.path.insert(0,'/tmp/gw'); import mr
from scipy import stats
from scipy.optimize import minimize_scalar
D=mr.load()
# add BMI and IGF-1 exposures from stage2 extraction
ref=pd.read_csv('/tmp/tools/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str}).set_index('rsid')[['a1','a2']]
def refck(d):
    d=d.join(ref,on='SNP'); d=d[((d.ea==d.a1)&(d.oa==d.a2))|((d.ea==d.a2)&(d.oa==d.a1))]; return d.drop(columns=['a1','a2'])
x=pd.read_csv('/tmp/gw/p5/bmi.tsv',sep='\t'); D['bmi']=mr.std(x,'A1','A2','Freq1.Hapmap','b','se','p'); D['bmi']['n']=233572
x=pd.read_csv('/tmp/gw/p5/igf.tsv',sep='\t'); x.columns=[c.strip() for c in x.columns]
D['igf']=refck(mr.std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')); D['igf']['n']=350000
clump={'biot':'/tmp/gw/biot_clump.clumped','psa':'/tmp/gw/psa_clump.clumped','bmi':'/tmp/gw/p2/bmi_clump.clumped','igf':'/tmp/gw/p2/igf_clump.clumped'}
def joint(names):
    u=pd.concat([pd.read_csv(clump[t],sep=r'\s+')[['SNP','P']] for t in names]).groupby('SNP',as_index=False).P.min()
    u.to_csv('/tmp/gw/p5/j.txt',sep='\t',index=False)
    subprocess.run(['/tmp/tools/plink','--bfile','/tmp/tools/EUR','--clump','/tmp/gw/p5/j.txt','--clump-p1','5e-8','--clump-p2','1','--clump-r2','0.001','--clump-kb','10000','--out','/tmp/gw/p5/j','--silent','--memory','2000'],check=True)
    return set(pd.read_csv('/tmp/gw/p5/j.clumped',sep=r'\s+').SNP)
def mvmr(m):
    X=np.column_stack([m.bx0,m.bx1]); w=1/m.sy.values**2; XtW=X.T*w
    b=np.linalg.solve(XtW@X,XtW@m.by.values); res=m.by.values-X@b; k=len(m)
    phi=max(1,np.sum(w*res**2)/(k-2)); se=np.sqrt(np.diag(np.linalg.inv(XtW@X))*phi)
    return b,se,k
rows=[]
for exp,lab in [('biot','Bioavailable testosterone (raises PSA)'),('bmi','BMI (lowers PSA)'),('igf','IGF-1 (no PSA effect)')]:
    inst=set(pd.read_csv(clump[exp],sep=r'\s+').SNP)
    uni=mr.harmonise(D[exp][D[exp].SNP.isin(inst)],D['prac'])
    b0,s0,_,_=mr.ivw(uni.b_x.values,uni.b_y.values,uni.se_y.values)
    snps=joint([exp,'psa']); base=D[exp][D[exp].SNP.isin(snps)]
    m=mr.harmonise(base,D['prac']).rename(columns={'b_x':'bx0','se_x':'sx0','b_y':'by','se_y':'sy'})[['SNP','bx0','sx0','by','sy']]
    h=mr.harmonise(base,D['psa']).rename(columns={'b_y':'bx1','se_y':'sx1'})[['SNP','bx1','sx1']]
    m=m.merge(h,on='SNP'); b,se,k=mvmr(m)
    rows.append(dict(exposure=lab,uni_OR=np.exp(b0),uni_p=2*stats.norm.sf(abs(b0/s0)),
                     adj_OR=np.exp(b[0]),adj_lo=np.exp(b[0]-1.96*se[0]),adj_hi=np.exp(b[0]+1.96*se[0]),
                     adj_p=2*stats.norm.sf(abs(b[0]/se[0])),change_pct=100*(b[0]-b0)/abs(b0),nsnp=k))
R=pd.DataFrame(rows); R.to_csv('/tmp/gw/p5/s2_collider_negative_control.tsv',sep='\t',index=False)
pd.set_option('display.width',200); print(R.round(3).to_string(index=False))
