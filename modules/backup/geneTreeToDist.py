from ete3 import Tree
import click, numpy as np, tempfile, os
from multiprocessing import Pool

def getDist(data) :
    names, fnames = data
    dist = np.zeros([len(names), len(names), 2], dtype=float)
    for fn in fnames :
        tre = Tree(fn, format=1)
        for node in tre.traverse('postorder') :
            if node.is_leaf() :
                assert node.name in names, 'tip name "{0}" is not present in the tip list'.format(node.name)
                node.d = {names[node.name]:0}
            else :
                for i, c1 in enumerate(node.children) :
                    for c2 in node.children[:i] :
                        for d1, l1 in c1.d.items() :
                            for d2, l2 in c2.d.items() :
                                dist[d1, d2, 0] += (l1+l2+1.)/len(names)
                                dist[d1, d2, 1] += 1.
                x = np.log(len(node.children))/np.log(2.)
                node.d = { dd:lvl+x for c in node.children for dd, lvl in c.d.items() }
    dist[:, :, 0] += dist[:, :, 0].T
    dist[:, :, 1] += dist[:, :, 1].T
    npz_file = tempfile.NamedTemporaryFile(dir='.', suffix='.npz', delete=False)
    npz_name = npz_file.name
    npz_file.close()
    np.savez_compressed(npz_name, dist=dist)
    return npz_name

@click.command()
@click.option('-t', '--tips', required=True)
@click.option('-o', '--outfile', required=True)
@click.argument('genetrees', nargs=-1)
def main(tips, outfile, genetrees) :
    pool = Pool(20)
    names = {}
    with open(tips) as fin :
        for line in fin :
            tip = line.strip().split()[0]
            if len(tip) :
                names[tip] = len(names)
    distances = np.zeros([len(names), len(names), 2], dtype=float)
    for npz_file in pool.imap_unordered(getDist, [ [names, genetrees[ite::20]] for ite in np.arange(20) ]) :
        dist = np.load(npz_file)['dist']
        os.unlink(npz_file)
        distances += dist

    distances[:, :, 0] += 0.0005
    distances[:, :, 1] += 0.0010
    distances = distances[:, :, 0] / distances[:, :, 1]
    np.fill_diagonal(distances, 0.)
    
    with open(outfile, 'wt') as fout :
        fout.write('    {0}\n'.format(distances.shape[0]))
        for (n, _), dist in zip(sorted(names.items(), key=lambda n:n[1]), distances) :
            fout.write('{0} {1}\n'.format(n, ' '.join([ '{:.12f}'.format(d) for d in dist.tolist()] )))
    return

if __name__ == '__main__' :
    main()