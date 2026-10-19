import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
import pandas as pd,gzip,sys
ids=set(pd.read_csv(f'{GW}/p5/union.txt',sep='\t').SNP)
bim=pd.read_csv(f'{TOOLS}/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str})
b=bim[bim.rsid.isin(ids)]; keys=dict(zip(b.chr+':'+b.pos.astype(str),b.rsid))
f,out,mode=sys.argv[1:4]
with gzip.open(f,'rt') as fh, open(out,'w') as o:
    h=fh.readline().rstrip('\r\n'); H=h.split('\t'); o.write(h+'\tSNP\n')
    for l in fh:
        x=l.rstrip('\r\n').split('\t')
        if mode=='rs0': s=x[0] if x[0] in ids else None
        elif mode=='key': s=keys.get(x[H.index('chromosome')]+':'+x[H.index('base_pair_location')])
        elif mode=='igf': s=x[H.index('rsid')] if x[H.index('rsid')] in ids else (x[H.index('variant_id')] if x[H.index('variant_id')] in ids else None)
        if s: o.write('\t'.join(x)+'\t'+s+'\n')
