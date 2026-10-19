import os as _os
GW = _os.environ.get('GW_WORK', '/tmp/gw')
TOOLS = _os.environ.get('GW_TOOLS', '/tmp/tools')
import gzip,sys,pandas as pd
bim=pd.read_csv(f'{TOOLS}/EUR.bim',sep='\t',header=None,names=['chr','rsid','cm','pos','a1','a2'],dtype={'chr':str})
bim=bim[bim.a1.str.len().eq(1)&bim.a2.str.len().eq(1)]
rs=set(bim.rsid); key=dict(zip(bim.chr+':'+bim.pos.astype(str),bim.rsid))
def run(f,out,pcol,mode):
    n=0
    with gzip.open(f,'rt') as fh, open(out,'w') as o:
        H=fh.readline().rstrip('\r\n').split('\t'); o.write('SNP\tP\n'); pi=H.index(pcol)
        for l in fh:
            x=l.rstrip('\r\n').split('\t')
            try: p=float(x[pi])
            except: continue
            if p>=5e-8: continue
            if mode=='bmi': s=x[0] if x[0] in rs else None
            elif mode=='igf': s=x[H.index('rsid')] if x[H.index('rsid')] in rs else (x[H.index('variant_id')] if x[H.index('variant_id')] in rs else None)
            else: s=key.get(x[0]+':'+x[1])
            if s: o.write(f'{s}\t{p}\n'); n+=1
    print(out,n)
job=sys.argv[1]
if job=='bmi': run('bmi.gz','bmi_clumpin.txt','p','bmi')
if job=='igf': run('igf1.h.tsv.gz','igf_clumpin.txt','p_value','igf')
if job=='fi': run('fi.tsv.gz','fi_clumpin.txt','p_value','fi')
