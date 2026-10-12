import gzip,sys
f,out,chrc,posc=sys.argv[1:5]; regs=[('12',2408330,3408330),('2',31482811,32482811)]
with gzip.open(f,'rt') as fh, open(out,'w') as o:
    h=fh.readline().rstrip('\r\n'); o.write(h+'\n'); H=h.split('\t'); ci,pi=H.index(chrc),H.index(posc)
    for l in fh:
        x=l.rstrip('\r\n').split('\t')
        for c,a,b in regs:
            if x[ci]==c and a<=int(x[pi])<=b: o.write('\t'.join(x)+'\n'); break
