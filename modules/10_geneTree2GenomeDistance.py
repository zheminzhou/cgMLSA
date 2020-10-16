
import os, sys, click, numpy as np, re, tempfile
from numba import njit
from ete3 import Tree
from multiprocessing import Pool
try :
    from .configure import externals, uopen
except :
    from configure import externals, uopen

@njit('void(f4[:, :], f8[:,:], f8[:,:], f4[:,:])')
def assign(dist, dd1, dd2, bounds) :
    for l1, d1 in dd1 :
        for l2, d2 in dd2 :
            if l1 > l2 :
                idx = int(l1 * (l1 - 1) / 2 + l2)
            else :
                idx = int(l2 * (l2 - 1) / 2 + l1)
            if not bounds.size or (bounds[idx, 0] <= d1+d2+5e-5 <= bounds[idx, 1]) :
                dist[idx, 0] += d1+d2
                dist[idx, 1] += 1


def treeDist(data) :
    treefiles, leaves, bounds = data
    dist = np.zeros([int(len(leaves) * (len(leaves) - 1) / 2), 2], dtype=np.float32)
    if bounds == None :
        bounds = np.zeros([0, 2], dtype=np.float32)
    else :
        bounds = np.load(bounds)['bounds']

    for treefile in treefiles :
        tre = Tree(treefile)
        for node in tre.traverse('postorder') :
            if node.is_leaf() :
                node.d = {leaves[node.name]:node.dist}
            else :
                for i1, c1 in enumerate(node.children) :
                    for c2 in node.children[:i1] :
                        assign(dist, np.array(list(c1.d.items())), np.array(list(c2.d.items())), bounds)
                node.d = { l:d+node.dist for c in node.children for l, d in c.d.items() }
    npz_file = tempfile.NamedTemporaryFile(dir='.', suffix='.npz', delete=False)
    npz_name = npz_file.name
    npz_file.close()
    np.savez_compressed(npz_name, dist=dist)
    return npz_name


@njit('void(f4[:, :], f8[:,:], f8[:,:])')
def assign_bin(dist, dd1, dd2) :
    for l1, d1 in dd1 :
        for l2, d2 in dd2 :
            if l1 > l2 :
                idx = int(l1 * (l1 - 1) / 2 + l2)
            else :
                idx = int(l2 * (l2 - 1) / 2 + l1)
            v = int(-6*np.log2(d1+d2+5e-4)+14)
            if v < 0 :
                v = 0
            elif v > 79 :
                v = 79
            dist[idx, v] += 1


def treeDistBin(data) :
    treefiles, leaves = data
    dist = np.zeros([int(len(leaves) * (len(leaves) - 1) / 2), 80], dtype=np.float32)

    for treefile in treefiles :
        tre = Tree(treefile)
        for node in tre.traverse('postorder') :
            if node.is_leaf() :
                node.d = {leaves[node.name]:node.dist}
            else :
                for i1, c1 in enumerate(node.children) :
                    for c2 in node.children[:i1] :
                        assign_bin(dist, np.array(list(c1.d.items())), np.array(list(c2.d.items())))
                node.d = { l:d+node.dist for c in node.children for l, d in c.d.items() }
    npz_file = tempfile.NamedTemporaryFile(dir='.', suffix='.npz', delete=False)
    npz_name = npz_file.name
    npz_file.close()
    np.savez_compressed(npz_name, dist=dist)
    return npz_name

@click.command()
@click.argument('trees', nargs=-1)
@click.option('-o', '--outfile', help='output file')
@click.option('-O', '--outlier', help='exclude outliers', default=False, is_flag=True)
def main(trees, outfile, outlier) :
    pool = Pool(10)
    leaves = set([])
    for ti, tree in enumerate(trees) :
        leaves |= set(re.findall('[(,]([^(),]+):', open(tree).read()))
    leaves = { n:i for i, n in enumerate(sorted(leaves))}

    if outlier :
        distance = np.zeros([int(len(leaves)*(len(leaves)-1)/2), 80], dtype=np.float32)
        for dist_npz in pool.imap_unordered(treeDistBin, [ (trees[i::20], leaves) for i in np.arange(20) ]) :
            dist = np.load(dist_npz)['dist']
            os.unlink(dist_npz)
            distance += dist
        distance = np.cumsum(distance, 1)
        distance = distance/distance[:, -1:]
        quantiles = -np.vstack([ np.argmax(distance >= 0.75, 1),  np.argmax(distance >= 0.25, 1)]).astype(np.float32)+13.5
        idx = (quantiles[0] == quantiles[1])
        quantiles[:, idx] += [[-0.5], [0.5]]
        quantiles[0, quantiles[0] <= -80 + 0.6] = -300
        quantiles = np.power(2, quantiles/6.)
        dqi = quantiles[1] - quantiles[0]
        quantiles += [-dqi*1.5, dqi*1.5]
        np.savez_compressed(outfile+'.data.npz', bounds=quantiles.T)
        quantiles = outfile+'.data.npz'
    else :
        quantiles = None

    distance = np.zeros([int(len(leaves)*(len(leaves)-1)/2), 2], dtype=np.float32)
    for dist_npz in pool.imap_unordered(treeDist, [ (trees[i::20], leaves, quantiles) for i in np.arange(20) ]) :
        dist = np.load(dist_npz)['dist']
        os.unlink(dist_npz)
        distance += dist
    distance = distance[:, 0]/distance[:, 1]
    if outlier:
        bounds = np.load(outfile+'.data.npz')['bounds']
        np.savez_compressed(outfile + '.data.npz', \
                            names=np.array([n for n, i in sorted(leaves.items(), key=lambda x:x[1])]), \
                            bounds=bounds, dist=distance)
        del bounds
    else :
        np.savez_compressed(outfile + '.data.npz', \
                            names=np.array([n for n, i in sorted(leaves.items(), key=lambda x:x[1])]), \
                            dist=distance)

    with uopen(outfile, 'w') as fout :
        distances2 = np.zeros([len(leaves), len(leaves)], dtype=np.float32)
        for i, n in enumerate(leaves) :
            dist = distance[ int(i*(i-1)/2):int(i*(i+1)/2) ]
            distances2[i, :i] = dist
            distances2[:i, i] = dist
        fout.write('    {0}\n'.format(len(leaves)))
        for n, i in sorted(leaves.items(), key=lambda x:x[1]) :
            fout.write('{0} {1}\n'.format(n, ' '.join([ '{:.8f}'.format(d) for d in distances2[i].tolist()] )))


if __name__ == '__main__' :
  main()
