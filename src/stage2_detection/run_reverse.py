import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
import pandas as pd,numpy as np,sys
sys.path.insert(0,GW); import mr
bim=pd.read_csv(f'{GW}/union_snps.tsv',sep='\t')  # not needed
ref=pd.read_csv(f'{TOOLS}/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str})
x=pd.read_csv('p2/r_prac.tsv',sep='\t'); x.columns=[c.strip() for c in x.columns]
p=mr.std(x,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')
r=ref[ref.rsid.isin(p.SNP)].set_index('rsid')[['a1','a2']]
p=p.join(r,on='SNP'); p=p[((p.ea==p.a1)&(p.oa==p.a2))|((p.ea==p.a2)&(p.oa==p.a1))].drop(columns=['a1','a2']); p['n']=140254
R=[];M=[];H=[]
for t,n in [('biot',178782),('tt',194453),('shbg',180726)]:
    y=pd.read_csv(f'p2/r_{t}.tsv',sep='\t'); o=mr.std(y,'effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value'); o['n']=n
    lg={}; m=mr.harmonise(p,o,lg); lg['outcome']=t; H.append(lg)
    res,meta=mr.run(m,f'prostate cancer liability->{t}'); R+=res; M.append(meta)
    tx=m.b_x/m.se_x; ty=m.b_y/m.se_y; rx=tx**2/(tx**2+m.n_x-2); ry=ty**2/(ty**2+m.n_y-2); st=(rx>ry).values
    res,meta=mr.run(m[st],f'prostate cancer liability->{t} | Steiger'); meta['steiger_removed']=int((~st).sum()); R+=res; M.append(meta)
pd.DataFrame(R).to_csv('p2/reverse_results.tsv',sep='\t',index=False); pd.DataFrame(M).to_csv('p2/reverse_diag.tsv',sep='\t',index=False)
pd.set_option('display.width',220)
r=pd.DataFrame(R); r['beta_per_doubling']=r.beta*np.log(2)
print(pd.DataFrame(H).to_string(index=False)); print(r.round(4).to_string(index=False)); print(pd.DataFrame(M).round(4).to_string(index=False))
