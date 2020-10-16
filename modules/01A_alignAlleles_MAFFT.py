import os, sys, click, numpy as np, pandas as pd, subprocess
from configure import readFasta, uopen, rc, externals
from multiprocessing import Pool

mafft = externals['mafft']
def alignAGene(data) :
    fin, fout = data
    p = subprocess.Popen('{0} --thread 1 --auto {1} | gzip > {2}'.format(mafft, fin, fout), stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()
    

def processGenes(data, outdir, cnt) :
    if len(data) < cnt :
        return data
    inputs = []
    for g, alleles in data.items() :
        fin = os.path.join(outdir, '{0}.alleles.fasta'.format(g))
        fout = os.path.join(outdir, '{0}.aln.fasta.gz'.format(g))
        num = 0
        with open(fin, 'wt') as i :
            for n, s in alleles.items() :
                s = '\n'.join(s)
                num += len(s)
                i.write('>{0}\n{1}\n'.format(n, s))
        inputs.append([num*len(alleles), fin, fout])
    list(pool.map(alignAGene, [ d[1:] for d in sorted(inputs, reverse=True) ]))
    for d in inputs :
        os.unlink(d[1])
        
    return {}

@click.command()
@click.option('-i', '--input', help='allelic sequences (can in gzip format)')
@click.option('-o', '--outdir', help='folder storing output alignment')
@click.option('-p', '--n_proc', help='n_proc', default=10, type=int)
def alignAlleles(input, outdir, n_proc) :
    global pool
    pool = Pool(n_proc)
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)
    with uopen(input, 'r') as fin :
        data = {}
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split(' ', 1)[0]
                gene, allele = name.rsplit('_', 1)
                if gene not in data :
                    data = processGenes(data, outdir, 40)
                    data[gene] = {}
                data[gene][name] = []
            else :
                data[gene][name].append(line.strip())
    processGenes(data, outdir, 1)
    

pool = None
if __name__ == '__main__' :
    alignAlleles()
