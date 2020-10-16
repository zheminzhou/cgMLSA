import os, sys, click, numpy as np, pandas as pd, tempfile
from numba import njit
from collections import defaultdict, OrderedDict
from multiprocessing import Pool
try :
    from .configure import readFasta, uopen
except :
    from configure import readFasta, uopen

@njit('void(f4[:, :, :], i8[:], i8[:], f4[:])')
def range_assign(mat, x, y, v) :
    for i in x :
        for j in y:
            if i != j :
                mat[i, j] = v

def geneDist(data) :
    data, norm_file, ite, min_size = data
    normalize = np.load(norm_file)
    group_ids, normalize, genomes, groups, group_dist = normalize['group_ids'], \
                                                        normalize['normalize'], \
                                                        normalize['genomes'], \
                                                        normalize['groups'], \
                                                        normalize['group_dist']
    group_contents = [ (group_ids == idx) for idx in np.arange(np.max(group_ids)+1) ]

    outfile = '{0}.{1}.gz'.format(norm_file, ite)
    with uopen(outfile, 'wt') as fout :
        for d_i, (alleles, alignment) in enumerate(data) :
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
                dist[i, i, 1] = np.sum(p)

            distance = np.zeros([len(group_ids), len(group_ids), 2], dtype=np.float32)
            for i, ci in enumerate(convs):
                range_assign(distance, ci, ci, dist[i, i])
                for j, cj in enumerate(convs[:i]):
                    range_assign(distance, ci, cj, dist[i, j])
                    range_assign(distance, cj, ci, dist[i, j])

            for id, grp in enumerate(group_ids) :
                if grp >= 0 and distance[id, :, 1].sum() > 0 :
                    gdist = []
                    for grp2, ct in enumerate(group_contents) :
                        d = distance[id, ct]
                        d = d[d.T[1]>0]
                        v = np.sum(d, 0)
                        if d.shape[0] < min_size :
                            x1 = np.sum(pres, 1).mean() * (min_size - d.shape[0])/2
                            x2 = group_dist[grp, grp2] * x1
                            v += [x2, x1]
                        gdist.append(v[0]/v[1]+5e-5)
                    gdist = np.array(gdist)/normalize[id]
                    fout.write('{0}\t{1}\t|\t{2}\t|\t{3}\n'.format(alignment, genomes[id], \
                                                                   groups[grp], \
                                                                   '\t'.join([ '{:.8f}'.format(v) for v in gdist.tolist()])))
    return outfile

def prepare_norm(distfile, grpfile, others, min_size, outfile) :
    with uopen(distfile) as fin :
        n_genome = int(fin.readline().strip().split()[0])
        genomes = []
        glob_dist = np.zeros( [n_genome, n_genome], dtype=np.float32 )
        for id, line in enumerate(fin) :
            part = line.strip().split()
            genomes.append(part[0])
            glob_dist[id] = part[1:]
    glob_dist += 5e-5
    np.fill_diagonal(glob_dist, 0.)
    genome_ids = { g:i for i, g in enumerate(genomes) }
    # assign groups
    tmp_rec = {}
    with open(grpfile) as fin :
        for line in fin :
            part = line.strip().split()
            if part[0] in genome_ids :
                tmp_rec[part[0]] = part[1]
    groups = np.unique(list(tmp_rec.values()), return_counts=True)
    groups = [grp for grp, cnt in zip(groups[0].tolist(), groups[1].tolist()) if cnt >= min_size]
    groups = { grp:id for id, grp in enumerate(sorted(groups)) }
    group_ids = np.array([-1]*n_genome, dtype=int)
    for genome, grp in tmp_rec.items() :
        group_ids[genome_ids[genome]] = groups.get(grp, -1)

    if others and np.sum(group_ids == -1) >= min_size :
        groups['others'] = len(groups)
        group_ids[group_ids == -1] = groups['others']

    group_contents = [ (group_ids == idx) for idx in np.arange(np.max(group_ids)+1) ]

    group_dist = np.zeros([len(groups), len(groups)], dtype=np.float32)
    for i, gc1 in enumerate(group_contents) :
        v = np.sum(glob_dist[gc1][:, gc1])/gc1.sum()/(gc1.sum()-1)
        group_dist[i, i] = v
        for j, gc2 in enumerate(group_contents[:i]) :
            v = np.mean(glob_dist[gc1][:, gc2])
            group_dist[j, i] = v
            group_dist[i, j] = v

    normalize = np.zeros([len(genomes), len(groups)])
    for i, dist in enumerate(glob_dist) :
        if group_ids[i] != -1 :
            jd = group_ids[i]
            normalize[i] = np.array([dist[gc].mean() if id != jd else dist[gc].sum()/(gc.sum()-1) \
                                     for id, gc in enumerate(group_contents)])/group_dist[group_ids[i]]

    np.clip(normalize, 0.01, 100., out=normalize)
    groups = np.array([g for g, i in sorted(groups.items(), key=lambda g: g[0])])
    np.savez_compressed(outfile+'.normalize.npz', \
                        genomes=np.array(genomes), \
                        group_ids=group_ids, \
                        normalize=normalize, \
                        groups=groups, \
                        group_dist=group_dist)
    return genome_ids, groups, outfile+'.normalize.npz'


@click.command()
@click.option('-p', '--profile', help='profile', required=True)
@click.option('-g', '--grpfile', help='group', required=True)
@click.option('-d', '--distfile', help='distance file', required=True)
@click.option('-o', '--outfile', help='outfile', required=True)
@click.option('-m', '--min_size', help='minimum size of a group', default=5, type=int)
@click.option('-O', '--others', help='label to assign ungrouped genomes to "others"', is_flag=True, default=False)
@click.argument('alignments', nargs=-1)
def main(profile, outfile, grpfile, distfile, min_size, others, alignments) :
    pool = Pool(10)
    # get genome names & distances
    genome_ids, groups, norm_file = prepare_norm(distfile, grpfile, others, min_size, outfile)

    data = pd.read_csv(profile, delimiter='\t', header=0, dtype=str)
    loci = data.columns
    data = data.values
    idx = np.argsort([genome_ids[g] for g in data.T[0]])
    data = data[idx]

    locus_files = []
    for fn in alignments :
        with uopen(fn) as fin :
            locus_files.append(fin.readline()[1:].strip().split()[0].rsplit('_', 1)[0])
    inputs = [ [allele, align] for allele, align in zip([data.T[loci == locus][0] for locus in locus_files], alignments)]
    stepwise = min(20, len(inputs))

    with uopen(outfile, 'wt') as fout :
        for grp in groups :
            fout.write('@GRP {0}\n'.format(grp))
        for dist_gz in pool.imap_unordered(geneDist, [ [inputs[ite::stepwise], norm_file, ite, min_size] \
                                                       for ite in np.arange(stepwise) ]) :
            with uopen(dist_gz) as fin :
                for line in fin :
                    fout.write(line)
            os.unlink(dist_gz)

if __name__ == '__main__' :
    main()
