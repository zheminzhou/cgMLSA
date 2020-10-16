import click, numpy as np, sys, gzip, numba as nb

@nb.njit()
def getDist(mat, dist) :
    l1 = mat.shape[0]
    l2 = mat.shape[1]
    for i in range(1, l1) :
        mi = mat[i]
        for j in range(i) :
            mj = mat[j]
            d = 0
            for k in range(l2) :
                if mi[k] != mj[k] and mi[k] != 45 and mj[k] != 45 :
                    d += 1
            dist[i, j] = d
            dist[j, i] = d
    return

@click.command()
@click.option('-r', '--rtab', help='PEPPAN rtab')
@click.option('-f', '--fmt', help='output format. phylip [default] or pairwise', default='phylip')
def main(rtab, fmt) :
    data = []
    with gzip.open(rtab, 'rt') as fin :
        for lId, line in enumerate(fin) :
            if line.startswith('##') :
                continue
            elif line.startswith('#') or lId == 0 :
                names = line.strip().split('\t')
                name_cols = np.array([ not n.startswith('#') for n in names ])
                name_cols[0] = False
                names = np.array(names)[name_cols]
            else :
                part = line.strip().split('\t')
                p = np.array(part)[name_cols]
                p[p != '0'] = '1'
                if max([len(pp) for pp in p]) == 1 :
                    data.append([ord(pp) for pp in p])

    data = np.array(data).T
    dist = np.zeros([data.shape[0], data.shape[0]], dtype=int)
    getDist(data, dist)
    if fmt == 'phylip' :
        sys.stdout.write('    {0}\n'.format(dist.shape[0]))
        for n, dd in zip(names, dist) :
            sys.stdout.write('{0} {1}\n'.format(n, ' '.join([ '{:.2f}'.format(d) for d in dd.tolist()] )))
    else :
        for i, di in enumerate(dist) :
            for j, d in enumerate(di[:i]) :
                sys.stdout.write('{0}\t{1}\t{2}\n'.format(names[j], names[i], d))


if __name__ == '__main__' :
    main()
