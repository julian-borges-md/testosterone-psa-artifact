import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
import pandas as pd,gzip,sys
ids=set()
for t in ['biot','shbg','tt','psa']:
    c=pd.read_csv(f'{t}_clump.clumped',sep=r'\s+'); ids|=set(c.SNP)
bim=pd.read_csv(f'{TOOLS}/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str})
b=bim[bim.rsid.isin(ids)].copy(); b['key']=b.chr+':'+b.pos.astype(str)
b.to_csv('union_snps.tsv',sep='\t',index=False)
keys=dict(zip(b.key,b.rsid)); rs=set(b.rsid)
print('union',len(rs))
def stream(f,out,keyfun):
    with gzip.open(f,'rt') as fh, open(out,'w') as o:
        h=fh.readline().rstrip('\r\n'); o.write(h+'\tSNP\n'); H=h.split('\t')
        for l in fh:
            x=l.rstrip('\r\n').split('\t'); s=keyfun(H,x)
            if s: o.write('\t'.join(x)+'\t'+s+'\n')
def byrs(col):
    def f(H,x):
        v=x[H.index(col)]; return v if v in rs else None
    return f
def bykey(c,p):
    def f(H,x):
        return keys.get(x[H.index(c)]+':'+x[H.index(p)])
    return f
def fg(H,x):
    for v in x[H.index('rsids')].split(','):
        if v in rs: return v
job=sys.argv[1]
if job=='biot': stream('biot_m.tsv.gz','biot_u.tsv',byrs('variant_id'))
if job=='shbg': stream('shbg_m.tsv.gz','shbg_u.tsv',byrs('variant_id'))
if job=='tt': stream('tt_m.tsv.gz','tt_u.tsv',byrs('variant_id'))
if job=='psa': stream('psa.tsv.gz','psa_u.tsv',bykey('chromosome','base_pair_location'))
if job=='prac': stream('pca_practical.tsv.gz','prac_u.tsv',bykey('chromosome','base_pair_location'))
if job=='fg': stream('pca_finngen.gz','fg_u.tsv',fg)
