import numpy as np, sys, subprocess, urllib, json, pandas as pd, gzip

if __name__ == '__main__' :
    seq = {}
    for fname in sys.argv[2:] :
        with gzip.open(fname) as fin :
            for line in fin :
                if line[0] == '>' :
                    name = line[1:].strip().split()[0]
                    seq[name] = []
                else :
                    seq[name].append( line.strip() )
        zero_name = name.rsplit('_', 1)[0] + '_0'
        seq[zero_name] = ['-'] * sum([len(s) for s in seq[name]])
    for n, s in seq.iteritems() :
        seq[n] = ''.join(s)
    data = pd.read_csv(sys.argv[1], delimiter='\t', header=None, dtype=str).as_matrix()
    for d in data[1:] :
        name = d[0]
        s = []
        for locus, allele in zip(data[0][2:], d[2:]) :
            key = '{0}_{1}'.format(locus, allele)
            if key in seq :
                s.append(seq[key])
            else :
                s.append(seq[locus+'_0'])
        print '>{0}\n{1}'.format(name, ''.join(s))
