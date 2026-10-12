import pandas as pd,numpy as np,sys,subprocess
sys.path.insert(0,'/tmp/gw'); import mr
from scipy import stats
ref=pd.read_csv('/tmp/tools/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str}).set_index('rsid')[['a1','a2']]
def refcheck(d):
    d=d.join(ref,on='SNP'); d=d[((d.ea==d.a1)&(d.oa==d.a2))|((d.ea==d.a2)&(d.oa==d.a1))]; return d.drop(columns=['a1','a2'])
L={}
x=pd.read_csv('mv_biot.tsv',sep='\t'); L['biot']=mr.std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value'); L['biot']['n']=178782
x=pd.read_csv('mv_bmi.tsv',sep='\t'); L['bmi']=mr.std(x,'A1','A2','Freq1.Hapmap','b','se','p'); L['bmi']['n']=233572
x=pd.read_csv('mv_igf.tsv',sep='\t'); x.columns=[c.strip() for c in x.columns]; L['igf']=refcheck(mr.std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')); L['igf']['n']=350000
x=pd.read_csv('mv_fi.tsv',sep='\t'); L['fi']=refcheck(mr.std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')); L['fi']['n']=151013
x=pd.read_csv('mv_prac.tsv',sep='\t'); x.columns=[c.strip() for c in x.columns]; L['prac']=refcheck(mr.std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')); L['prac']['n']=140254
x=pd.read_csv('mv_fg.tsv',sep='\t'); L['fg']=mr.std(x,'alt','ref','af_alt','beta','sebeta','pval'); L['fg']['n']=np.nan
clump={'biot':'/tmp/gw/biot_clump.clumped','bmi':'bmi_clump.clumped','igf':'igf_clump.clumped','fi':'fi_clump.clumped'}
def joint(names):
    u=pd.concat([pd.read_csv(clump[t],sep=r'\s+')[['SNP','P']] for t in names]).groupby('SNP',as_index=False).P.min()
    u.to_csv('jin.txt',sep='\t',index=False)
    subprocess.run(['/tmp/tools/plink','--bfile','/tmp/tools/EUR','--clump','jin.txt','--clump-p1','5e-8','--clump-p2','1','--clump-r2','0.001','--clump-kb','10000','--out','jj','--silent','--memory','2000'],check=True)
    return set(pd.read_csv('jj.clumped',sep=r'\s+').SNP)
def build(names,out):
    snps=joint(names); base=L[names[0]][L[names[0]].SNP.isin(snps)]
    m=mr.harmonise(base,L[out]).rename(columns={'b_x':'bx0','se_x':'sx0','b_y':'by','se_y':'sy'})[['SNP','bx0','sx0','by','sy']]
    for i,t in enumerate(names[1:],1):
        h=mr.harmonise(base,L[t]).rename(columns={'b_y':f'bx{i}','se_y':f'sx{i}'})[['SNP',f'bx{i}',f'sx{i}']]
        m=m.merge(h,on='SNP')
    return m,len(snps)
def mvmr(m,k_exp):
    X=np.column_stack([m[f'bx{i}'] for i in range(k_exp)]); w=1/m.sy.values**2; XtW=X.T*w
    b=np.linalg.solve(XtW@X,XtW@m.by.values); res=m.by.values-X@b; k=len(m)
    phi=max(1,np.sum(w*res**2)/(k-k_exp)); se=np.sqrt(np.diag(np.linalg.inv(XtW@X))*phi)
    cf=[]
    for j in range(k_exp):
        oth=[i for i in range(k_exp) if i!=j]; Xo=np.column_stack([m[f'bx{i}'] for i in oth]); wj=1/m[f'sx{j}'].values**2
        d=np.linalg.lstsq(Xo*np.sqrt(wj)[:,None],m[f'bx{j}'].values*np.sqrt(wj),rcond=None)[0]
        resj=m[f'bx{j}'].values-Xo@d
        var=m[f'sx{j}'].values**2+sum((d[ii]**2)*m[f'sx{oth[ii]}'].values**2 for ii in range(len(oth)))
        cf.append(np.sum(resj**2/var)/(k-len(oth)))
    return b,se,cf,k,phi
rows=[]
for out in ['prac','fg']:
    for names in [['biot','bmi'],['biot','igf'],['biot','fi'],['biot','bmi','igf','fi']]:
        m,nj=build(names,out); b,se,cf,k,phi=mvmr(m,len(names))
        for i,n in enumerate(names):
            rows.append(dict(outcome=out,model='+'.join(names),exposure=n,nsnp=k,beta=b[i],se=se[i],OR=np.exp(b[i]),lo=np.exp(b[i]-1.96*se[i]),hi=np.exp(b[i]+1.96*se[i]),p=2*stats.norm.sf(abs(b[i]/se[i])),condF=cf[i]))
uni=[]
for t in ['bmi','igf','fi']:
    ex=L[t][L[t].SNP.isin(set(pd.read_csv(clump[t],sep=r'\s+').SNP))]
    for out in ['prac','fg']:
        m=mr.harmonise(ex,L[out]); r,meta=mr.run(m,f'{t}->{out}'); uni+=r
r=pd.DataFrame(rows); r.to_csv('mvmr_metabolic.tsv',sep='\t',index=False)
u=pd.DataFrame(uni); u.to_csv('uni_metabolic.tsv',sep='\t',index=False)
pd.set_option('display.width',220); print(r.round(4).to_string(index=False)); print(u[u.method.isin(['IVW','Weighted median'])].round(4).to_string(index=False))
