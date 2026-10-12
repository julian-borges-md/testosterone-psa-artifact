import pandas as pd,numpy as np,mr,json
D=mr.load(); bim=pd.read_csv('union_snps.tsv',sep='\t',dtype={'chr':str}).set_index('rsid')
inst={t:set(pd.read_csv(f'{t}_clump.clumped',sep=r'\s+').SNP) for t in ['biot','shbg','tt','psa']}
def region(s,chr_,lo,hi):
    r=bim.loc[s]; return (r.chr==chr_)&(r.pos>=lo)&(r.pos<=hi)
def steiger(m):
    tx=m.b_x/m.se_x; ty=m.b_y/m.se_y
    nx=m.n_x; ny=m.n_y.fillna(m.n_y.median() if m.n_y.notna().any() else 300000)
    rx=tx**2/(tx**2+nx-2); ry=ty**2/(ty**2+ny-2); return rx>ry
R=[];M=[];H=[]
plan=[('biot','psa'),('shbg','psa'),('tt','psa'),('biot','prac'),('biot','fg'),('shbg','prac'),('tt','prac'),('psa','prac'),('psa','fg'),('psa','biot')]
for e,o in plan:
    ex=D[e][D[e].SNP.isin(inst[e])]
    lg={}; m=mr.harmonise(ex,D[o],lg); lg.update(exposure=e,outcome=o,instruments=len(ex)); H.append(lg)
    lab=f'{e}->{o}'
    r,mt=mr.run(m,lab); R+=r; M.append(mt)
    st=steiger(m); r,mt=mr.run(m[st.values],lab+' | Steiger'); R+=r; mt['steiger_removed']=int((~st).sum()); M.append(mt)
    excl=~(region(m.SNP,'19',50300000,52400000).values|region(m.SNP,'17',6500000,8600000).values)
    if excl.sum()<len(m):
        r,mt=mr.run(m[excl],lab+' | excl SHBG/KLK3 loci'); R+=r; mt['loci_removed']=int((~excl).sum()); M.append(mt)
pd.DataFrame(R).to_csv('out/mr_results.tsv',sep='\t',index=False)
pd.DataFrame(M).to_csv('out/mr_diagnostics.tsv',sep='\t',index=False)
pd.DataFrame(H).to_csv('out/harmonisation_log.tsv',sep='\t',index=False)
