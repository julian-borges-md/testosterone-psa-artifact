import pandas as pd,numpy as np
from scipy.special import logsumexp
a=pd.read_csv('biot_reg.tsv',sep='\t',dtype={'chromosome':str}); b=pd.read_csv('psa_reg.tsv',sep='\t',dtype={'chromosome':str})
for d in (a,b):
    d['effect_allele']=d.effect_allele.str.upper(); d['other_allele']=d.other_allele.str.upper()
    d.drop(d[(d.effect_allele.str.len()!=1)|(d.other_allele.str.len()!=1)].index,inplace=True)
    d['key']=d.chromosome+':'+d.base_pair_location.astype(str)
a=a[a.imputation_quality>=0.8]
m=a.merge(b,on='key',suffixes=('_a','_b'))
ok=((m.effect_allele_a==m.effect_allele_b)&(m.other_allele_a==m.other_allele_b))|((m.effect_allele_a==m.other_allele_b)&(m.other_allele_a==m.effect_allele_b))
m=m[ok].drop_duplicates('key')
def labf(beta,se,maf,N,W2=None):
    V=se**2
    if W2 is None:
        sdy2=np.nanmedian(V*2*N*maf*(1-maf)); W2=0.15**2*sdy2
    r=W2/(W2+V); z=beta/se
    return 0.5*(np.log(1-r)+r*z**2)
def coloc(l1,l2,p1=1e-4,p2=1e-4,p12=1e-5):
    lH1=np.log(p1)+logsumexp(l1); lH2=np.log(p2)+logsumexp(l2)
    l12=l1+l2; lH4=np.log(p12)+logsumexp(l12)
    lsum=logsumexp(l1)+logsumexp(l2)
    lH3=np.log(p1)+np.log(p2)+lsum+np.log1p(-np.exp(logsumexp(l12)-lsum))
    L=np.array([0,lH1,lH2,lH3,lH4]); return np.exp(L-logsumexp(L))
for c,lab in [('12','FKBP4 region chr12:2.41-3.41Mb'),('2','2p22 region chr2:31.48-32.48Mb')]:
    s=m[m.chromosome_a==c]
    maf=s.effect_allele_frequency_a.clip(0,1); maf=np.minimum(maf,1-maf); s=s[maf>0.005]; maf=maf[maf>0.005]
    l1=labf(s.beta_a.values,s.standard_error_a.values,maf.values,178782)
    l2=labf(s.beta_b.values,s.standard_error_b.values,maf.values,s.n.values.astype(float))
    pp=coloc(l1,l2)
    top=s.iloc[np.argmax(l1+l2)]
    print(lab,'nsnp',len(s),'PP H0-H4',np.round(pp,4),'lead shared',top.key,'p_biot %.1e p_psa %.1e'%(top.p_value_a,top.p_value_b))
