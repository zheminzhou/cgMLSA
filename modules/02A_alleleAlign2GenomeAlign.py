import numpy as np, sys, subprocess, urllib, json, pandas as pd, gzip, click

try:
    from .configure import readFasta, uopen
except:
    from configure import readFasta, uopen


@click.command()
@click.option('-p', '--profile')
@click.option('-s', '--subset', default=None)
@click.option('-o', '--output')
@click.argument('alleles', nargs=-1)
def main(output, profile, subset, alleles):
    '''
    alleles: aligned allelic sequences in multi-fasta format
    '''
    data = pd.read_csv(profile, delimiter='\t', header=None, dtype=str).values
    genomes = data.T[0, 1:]
    genes = data[0, 1:]
    data = data[1:, 1:]

    if subset:
        subsets = subset.split(',')
        index = np.in1d(genomes, subsets)
        genomes = genomes[index]
        data = data[index]
    #sequences = {n: [] for n in genomes}

    with uopen(output, 'w') as fout :
        for fname in alleles:
            seqs = readFasta(fname)
            n, s = list(seqs.items())[0]
            locus = n.rsplit('_', 1)[0]
            no_seq = '-' * len(s)
            loc_seqs = []
            if locus in genes:
                genome_alleles = data.T[genes == locus][0]
                for a in genome_alleles:
                    loc_seqs.append(list(seqs.get('{0}_{1}'.format(locus, a), no_seq)))
                loc_seqs = np.array(loc_seqs)
                loc_seqs = loc_seqs[:, np.any(loc_seqs != '-', 0)]
                for n, a, i in zip(genomes, loc_seqs,  genome_alleles) :
                    fout.write('>{0}:{1} {2}\n{3}\n'.format(n, locus, i, ''.join(a.tolist())))
                fout.write('=\n')

if __name__ == '__main__' :
    main()

