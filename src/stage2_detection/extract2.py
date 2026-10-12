import pandas as pd,gzip,sys
ids=set(pd.read_csv(sys.argv[1],sep=r'\s+').SNP)
bim=pd.read_csv('/tmp/tools/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str})
b=bim[bim.rsid.isin(ids)].copy(); b['key']=b.chr+':'+b.pos.astype(str); keys=dict(zip(b.key,b.rsid))
f,out,mode=sys.argv[2],sys.argv[3],sys.argv[4]
with gzip.open(f,'rt') as fh, open(out,'w') as o:
    h=fh.readline().rstrip('\r\n'); o.write(h+'\tSNP\n'); H=h.split('\t')
    for l in fh:
        x=l.rstrip('\r\n').split('\t')
        if mode=='rs':
            s=x[0] if x[0] in ids else None
        else:
            s=keys.get(x[H.index('chromosome')]+':'+x[H.index('base_pair_location')])
        if s: o.write('\t'.join(x)+'\t'+s+'\n')
