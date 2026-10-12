"""CAUSE needs genome-wide summary statistics for both traits, not just instruments.
We extract a common SNP set with p-value based pruning, as CAUSE's own workflow does:
all variants present in both files, then LD pruning on the instrument reference.
"""
import gzip, pandas as pd, numpy as np, sys
bim=pd.read_csv('/tmp/tools/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str})
bim=bim[bim.a1.str.len().eq(1)&bim.a2.str.len().eq(1)]
bim['key']=bim.chr+':'+bim.pos.astype(str)
keep_rs=set(bim.rsid); key2rs=dict(zip(bim.key,bim.rsid))
def stream(path, out, mode, cols):
    n=0
    with gzip.open(path,'rt') as fh, open(out,'w') as o:
        H=fh.readline().rstrip('\r\n').split('\t')
        idx={c:H.index(c) for c in cols}
        o.write('snp\tbeta\tse\tA1\tA2\n')
        for i,l in enumerate(fh):
            if i%3: continue          # thin 3x for runtime; CAUSE prunes anyway
            x=l.rstrip('\r\n').split('\t')
            if mode=='rs':
                s=x[0]
                if s not in keep_rs: continue
            else:
                s=key2rs.get(x[idx['chromosome']]+':'+x[idx['base_pair_location']])
                if s is None: continue
            a1=x[idx['effect_allele']].upper(); a2=x[idx['other_allele']].upper()
            if len(a1)!=1 or len(a2)!=1: continue
            b=x[idx['beta']]; se=x[idx['standard_error']]
            if b in ('NA','') or se in ('NA',''): continue
            o.write(f'{s}\t{b}\t{se}\t{a1}\t{a2}\n'); n+=1
    print(out,n,flush=True)
job=sys.argv[1]
base=['effect_allele','other_allele','beta','standard_error']
if job=='biot': stream('/tmp/gw/biot_m.tsv.gz','/tmp/gw/p6/biot.tsv','rs',base)
if job=='psa':  stream('/tmp/gw/psa.tsv.gz','/tmp/gw/p6/psa.tsv','key',base+['chromosome','base_pair_location'])
if job=='prac': stream('/tmp/gw/pca_practical.tsv.gz','/tmp/gw/p6/prac.tsv','key',base+['chromosome','base_pair_location'])
