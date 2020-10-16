import numpy as np, os, sys, subprocess, urllib, json, pandas as pd, gzip, click
from configure import uopen

@click.command()
@click.option('-p', '--profile', help='profile', required=True)
@click.option('-o', '--outdir', help='Folder for output fasta files. Default: current folder', default='.')
@click.argument('fasfiles', nargs=-1)
def main(profile, outdir, fasfiles) :
    data = pd.read_csv(profile, delimiter='\t', header=0, dtype=str)
    loci = data.columns
    data = data.values
    genomes = data.T[0, 1:]
    seqs = {}
    for i, (n, d) in enumerate(zip(loci[1:], data[1:, 1:].T)) :
        seqs[n] = {}
        for g, v in zip(genomes, d) :
            if v not in seqs[n] :
                seqs[n][v] = [g]
            else :
                seqs[n][v].append(g)
    fout = None
    for fname in fasfiles :
        with uopen(fname, 'r') as fin :
            toWrite = False
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    n, v = name.rsplit('_', 1)
                    n = n.split(':', 1)[-1]
                    if n in seqs and v in seqs[n] :
                        #print(n, v)
                        toWrite = True
                        outfile = os.path.join(outdir, '{0}.alleles.fas'.format(n))
                        if fout != None and fout.name != outfile :
                            fout.close()
                            fout = None
                        if fout == None :
                            fout = open(outfile, 'at')
                        fout.write('>{0}_{1} {2}\n'.format(n, v, ' '.join(seqs[n][v])))
                    else :
                        toWrite=False
                elif toWrite :
                    fout.write(line)
    if fout != None :
        fout.close()

if __name__ == '__main__' :
    main()