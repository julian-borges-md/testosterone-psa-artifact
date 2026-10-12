import pandas as pd
bim=pd.read_csv('/tmp/tools/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str})
bim=bim[bim.a1.str.len().eq(1)&bim.a2.str.len().eq(1)]
bim['key']=bim.chr+':'+bim.pos.astype(str)
log=[]
def qc(df,name,has_info):
    n0=len(df)
    df=df[df.effect_allele.str.len().eq(1)&df.other_allele.str.len().eq(1)]
    n1=len(df)
    maf=df.effect_allele_frequency.clip(upper=1); maf=maf.where(maf<=0.5,1-maf)
    df=df[maf>=0.01]; n2=len(df)
    if has_info: df=df[df.imputation_quality>=0.8]
    n3=len(df)
    log.append((name,n0,n1,n2,n3)); return df
for t in ['biot','shbg','tt']:
    d=pd.read_csv(f'{t}_sig.tsv',sep='\t',dtype={'chromosome':str})
    d=qc(d,t,True)
    d=d[d.variant_id.isin(set(bim.rsid))]
    log.append((t+'_inref',len(d)))
    d[['variant_id','p_value']].rename(columns={'variant_id':'SNP','p_value':'P'}).to_csv(f'{t}_clumpin.txt',sep='\t',index=False)
d=pd.read_csv('psa_sig.tsv',sep='\t',dtype={'chromosome':str})
d['effect_allele']=d.effect_allele.str.upper(); d['other_allele']=d.other_allele.str.upper()
d=qc(d,'psa',False)
d['key']=d.chromosome+':'+d.base_pair_location.astype(str)
m=d.merge(bim[['key','rsid','a1','a2']],on='key')
m=m[((m.effect_allele==m.a1)&(m.other_allele==m.a2))|((m.effect_allele==m.a2)&(m.other_allele==m.a1))]
m=m.drop_duplicates('key')
log.append(('psa_inref',len(m)))
m[['rsid','p_value']].rename(columns={'rsid':'SNP','p_value':'P'}).to_csv('psa_clumpin.txt',sep='\t',index=False)
for l in log: print(l)
