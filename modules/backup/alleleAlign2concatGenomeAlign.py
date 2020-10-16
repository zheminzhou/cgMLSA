import numpy as np, sys, subprocess, urllib, json, pandas as pd, gzip, click

try:
    from .configure import readFasta
except:
    from configure import readFasta


@click.command()
@click.option('-p', '--profile')
@click.option('-g', '--group', default=None)
@click.option('-o', '--output')
@click.argument('alleles', nargs=-1)
def main(output, profile, group, alleles):
    data = pd.read_csv(profile, delimiter='\t', header=None, dtype=str).values
    genomes = data.T[0, 1:]
    genes = data[0, 1:]
    data = data[1:, 1:]

    if group:
        groups = group.split(',')
        index = np.in1d(genomes, groups)
        genomes = genomes[index]
        data = data[index]
    sequences = {n: [] for n in genomes}

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
            loc_seqs = loc_seqs[:, np.all(loc_seqs != '-', 0)]
            for n, a in zip(genomes, loc_seqs) :
                sequences[n].append(''.join(a.tolist()))
    with gzip.open(output, 'wt') as fout:
        for n, s in sequences.items():
            s = ''.join(s)
            fout.write('>{0}\n{1}\n'.format(n, s))


if __name__ == '__main__' :
    main()