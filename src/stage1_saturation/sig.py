import gzip,sys
f,out,pcol=sys.argv[1],sys.argv[2],sys.argv[3]
with gzip.open(f,'rt') as fh, open(out,'w') as o:
    h=fh.readline().rstrip('\r\n').split('\t'); o.write('\t'.join(h)+'\n'); i=h.index(pcol)
    for l in fh:
        x=l.rstrip('\r\n').split('\t')
        try:
            if float(x[i])<5e-8: o.write('\t'.join(x)+'\n')
        except: pass
