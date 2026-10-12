import pandas as pd,numpy as np,mr
from scipy import stats
D=mr.load(); bim=pd.read_csv('union_snps.tsv',sep='\t',dtype={'chr':str}).set_index('rsid')
inst={t:set(pd.read_csv(f'{t}_clump.clumped',sep=r'\s+').SNP) for t in ['biot','shbg','tt','psa']}
# outliers in biot->psa
ex=D['biot'][D['biot'].SNP.isin(inst['biot'])]; m=mr.harmonise(ex,D['psa'])
b,s,Q,k=mr.ivw(m.b_x.values,m.b_y.values,m.se_y.values)
q=(m.b_x/m.se_y)**2*(m.b_y/m.b_x-b)**2; m['qp']=stats.chi2.sf(q,1); m['ratio']=m.b_y/m.b_x
m=m.join(bim[['chr','pos']],on='SNP')
o=m[m.qp<0.05/k].sort_values('qp')
print(o[['SNP','chr','pos','b_x','b_y','ratio','qp']].to_string(index=False))
m.to_csv('out/biot_psa_snps.tsv',sep='\t',index=False)
