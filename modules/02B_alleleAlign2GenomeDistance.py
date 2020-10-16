import os, sys, click, numpy as np, pandas as pd, tempfile
from numba import njit
from collections import defaultdict
from multiprocessing import Pool
try :
    from .configure import readFasta, uopen
except :
    from configure import readFasta, uopen

@njit('void(f4[:, :], i8[:], i8[:], f4[:], f4[:,:])')
def range_assign(mat, x, y, v, bounds) :
    for i in x :
        for j in y:
            if i != j :
                if i > j :
                    idx = int(i*(i-1)/2+j)
                else :
                    idx = int(j*(j-1)/2+i)
                if not bounds.size or ((v[1]+1.) * bounds[idx, 0] <= v[0]+0.5 <= (v[1]+1.) * bounds[idx, 1]) :
                    mat[idx, :] += v

@njit('void(f4[:, :], i8[:], i8[:], i4[:])')
def range_assign_bin(mat, x, y, v) :
    for i in x :
        for j in y:
            if i != j :
                if i > j :
                    idx = int(i * (i - 1)/2 + j)
                else :
                    idx = int(j * (j - 1)/2 + i)
                mat[idx, v[0]] += v[1]

def geneDist(d) :
    data, bounds, model = d
    if bounds :
        bounds = np.load(bounds)['bounds']
    else :
        bounds = np.zeros([0, 2], dtype=np.float32)
    alleles = data[0][0]
    distance = np.zeros([int(alleles.shape[0]*(alleles.shape[0]-1)/2.), 2], dtype=np.float32)
    for alleles, alignment in data :
        seqs = readFasta(alignment)
        names = np.array([ n.rsplit('_', 1)[-1] for n in seqs.keys()])

        tags = defaultdict(list)
        for j, v in enumerate(alleles):
            tags[v].append(j)
        convs = [np.array(tags.get(n, []), dtype=int) for n in names]

        seqs = np.array([ list(s) for s in seqs.values() ])
        pres = (seqs != '-')
        dist = np.zeros([len(names), len(names), 2], dtype=np.float32)
        for i, (s, p) in enumerate(zip(seqs, pres)) :
            iden = np.sum((seqs[:i] == s) & p, 1).astype(np.float32)
            shared = np.sum(pres[:i] & p, 1)
            dist[i, :i, 0] = shared - iden
            dist[i, :i, 1] = shared
            dist[i, i, 1] = np.sum(p)/2.
        dist[(dist[:, :, 1] < 0.01), :] = 0.01
        dist[(dist[:, :, 0] < 1e-8), 0] = 1e-8
        if model == 'jc' :
            raw_d = np.clip(dist[:, :, 0]/dist[:, :, 1], 1e-6, 0.745)
            raw_d = -3. / 4. * np.log(1 - 4 / 3. * raw_d)
            dist[:, :, 0] = dist[:, :, 1]*raw_d

        for i, ci in enumerate(convs):
            range_assign(distance, ci, ci, dist[i, i], bounds)
            for j, cj in enumerate(convs[:i]):
                range_assign(distance, ci, cj, dist[i, j], bounds)

    npz_file = tempfile.NamedTemporaryFile(dir='.', suffix='.npz', delete=False)
    npz_name = npz_file.name
    npz_file.close()
    np.savez_compressed(npz_name, dist=distance)
    return npz_name

def geneDistToBin(data) :
    data, bins, model = data
    alleles = data[0][0]
    distance = np.zeros([int(alleles.shape[0]*(alleles.shape[0]-1)/2), 90], dtype=np.float32)
    for alleles, alignment in data :
        seqs = readFasta(alignment)
        names = np.array([ n.rsplit('_', 1)[-1] for n in seqs.keys()])

        tags = defaultdict(list)
        for j, v in enumerate(alleles):
            tags[v].append(j)
        convs = [np.array(tags.get(n, []), dtype=int) for n in names]

        seqs = np.array([ list(s) for s in seqs.values() ])
        pres = (seqs != '-')
        dist = np.zeros([len(names), len(names), 2], dtype=np.float32)
        for i, (s, p) in enumerate(zip(seqs, pres)) :
            iden = np.sum((seqs[:i] == s) & p, 1).astype(np.int32)
            shared = np.sum(pres[:i] & p, 1)
            dist[i, :i, 0] = shared - iden
            dist[i, :i, 1] = shared
            dist[i, i, 1] = np.sum(p)
        dist[(dist[:, :, 1] < 0.01), :] = 0.01
        dist[(dist[:, :, 0] < 1e-8), 0] = 1e-8
        raw_d = dist[:, :, 0]/dist[:, :, 1]
        if model =='jc' :
            raw_d = np.clip(raw_d, 1e-6, 0.745)
            raw_d = -3./4.*np.log(1-4/3.*raw_d)
        dist[:, :, 0] = dist[:, :, 1]*raw_d + 0.5
        dist[:, :, 1] += 1.
        dist[:, :, 0] = -5*(np.log2(dist[:, :, 0]/dist[:, :, 1])-2)
        dist[:, :, 1] *=2
        np.fill_diagonal(dist[:, :, 1], np.diagonal(dist[:, :, 1])/2)

        dist = dist.astype(np.int32)
        np.clip(dist[:, :, 0], 0, 90-1, dist[:, :, 0])
        for i, ci in enumerate(convs):
            range_assign_bin(distance, ci, ci, dist[i, i])
            for j, cj in enumerate(convs[:i]):
                range_assign_bin(distance, ci, cj, dist[i, j])

    npz_file = tempfile.NamedTemporaryFile(dir='.', suffix='.npz', delete=False)
    npz_name = npz_file.name
    npz_file.close()
    np.savez_compressed(npz_name, dist=distance)
    return npz_name


@click.command()
@click.option('-p', '--profile', help='profile', required=True)
@click.option('-o', '--outfile', help='distance file', required=True)
@click.option('-O', '--outliers', help='remove outliers', is_flag=True, default=False)
@click.option('-m', '--model', help='raw or jc [default]', default='jc')
@click.argument('alignments', nargs=-1)
def main(profile, outfile, outliers, alignments, model) :
    pool = Pool(5)
    data = pd.read_csv(profile, delimiter='\t', header=0, dtype=str)
    loci = data.columns
    data = data.values
    genomes = data.T[0]

    locus_files = []
    for fn in alignments :
        with uopen(fn) as fin :
            locus_files.append(fin.readline()[1:].strip().split()[0].rsplit('_', 1)[0])
    
    inputs = list(zip([data.T[loci == locus][0] for locus in locus_files], alignments))

    if outliers :
        distances = np.zeros([int(genomes.size*(genomes.size-1)/2), 90], dtype=np.float32)
        bins = []
        for gId, dist_npz in enumerate(pool.imap_unordered(geneDistToBin, [ [inputs[ite::20], bins, model ] for ite in np.arange(20) ])) :
            dist = np.load(dist_npz)['dist']
            os.unlink(dist_npz)
            distances[:, :] += dist
        distances = np.cumsum(distances, 1)
        distances = distances/distances[:, -1:]
        quantiles = -np.vstack([ np.argmax(distances >= 0.75, 1),  np.argmax(distances >= 0.25, 1)]).astype(np.float32)-0.5
        idx = (quantiles[0] == quantiles[1])
        quantiles[:, idx] += [[-0.5], [0.5]]
        quantiles[0, quantiles[0] <= -90 + 0.6] = -300
        quantiles = np.power(2, quantiles/5.+2.)
        dqi = quantiles[1] - quantiles[0]
        quantiles += [-dqi*1.5, dqi*1.5]
        np.savez_compressed(outfile+'.data.npz', bounds=quantiles.T)
        quantiles = outfile+'.data.npz'
    else :
        quantiles = None
    distances = np.zeros([int(genomes.size*(genomes.size-1)/2), 2], dtype=np.float32)
    for gId, dist_npz in enumerate(pool.imap_unordered(geneDist, [ [inputs[ite::20], quantiles, model ] for ite in np.arange(20) ])) :
        dist = np.load(dist_npz)['dist']
        os.unlink(dist_npz)
        distances[:, :] += dist
    distances[:, 0] += 0.00005
    distances[:, 1] += 0.00010
    distances = distances[:, 0]/distances[:, 1]
    if outliers:
        bounds = np.load(outfile+'.data.npz')['bounds']
        np.savez_compressed(outfile + '.data.npz', names=genomes, bounds=bounds, dist=distances)
        del bounds
    else :
        np.savez_compressed(outfile + '.data.npz', names=genomes, dist=distances)

    with uopen(outfile, 'w') as fout :
        distances2 = np.zeros([genomes.size, genomes.size], dtype=np.float32)
        for i, n in enumerate(genomes) :
            dist = distances[ int(i*(i-1)/2):int(i*(i+1)/2) ]
            distances2[i, :i] = dist
            distances2[:i, i] = dist
        fout.write('    {0}\n'.format(len(genomes)))
        for i, n in enumerate(genomes) :
            fout.write('{0} {1}\n'.format(n, ' '.join([ '{:.8f}'.format(d) for d in distances2[i].tolist()] )))


if __name__ == '__main__' :
    main()
