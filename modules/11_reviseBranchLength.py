import sys, os, click, numpy as np, subprocess, re
from ete3 import Tree
try :
    from .configure import uopen, externals
except :
    from configure import uopen, externals

astral = externals['astral']
erable = externals['erable']

def assign_qs(tree) :
    tre = Tree(tree, format=1, quoted_node_names=True)
    tre.qs = np.array([1., 0., 0.])
    for node in tre.iter_descendants('postorder') :
        if not node.is_leaf() :
            qs = re.findall('\[q1=(.+);q2=(.+);q3=(.+)\]', node.name)
            try :
                node.qs = np.array([ float(q) for q in qs[0]])
                if np.isnan(node.qs[0]) :
                    node.qs = np.zeros(3, dtype=float)
            except Exception :
                node.qs = np.zeros(3, dtype=float)
    for node in tre.traverse('postorder') :
        if node.is_leaf() :
            node.d = (node.name,)
        else :
            node.d = tuple(sorted([ n for c in node.children for n in c.d ]))
    return tre

def assign_iqr(outfile, tree, genetree, n_split, n_thread) :
    if n_split <= 1 :
        geneTrees = [[1, genetree]]
    else :
        geneTrees = []
        with open(genetree) as fin :
            trees = fin.readlines()
        for ite in np.arange(n_split) :
            s_file = outfile + '.split.{0}'.format(ite)
            with open(s_file, 'wt') as fout :
                subtrees = trees[ite::n_split]
                fout.write('\n'.join(subtrees))
            geneTrees.append([len(subtrees), s_file])
    # run astral
    dqs_sum = {}
    tre = None
    for i, (num, tfile) in enumerate(geneTrees) :
        p = subprocess.Popen('{0} -q {1} -i {2} -t 8 -o {3} -T {4}'.format(astral, tree, tfile, outfile+'.astral', n_thread).split(), \
                             stdout=subprocess.PIPE)
        p.communicate()

        tre = assign_qs('{0}.astral'.format(outfile))
        if tfile != genetree :
            os.unlink(tfile)
        #os.unlink('{0}.astral'.format(outfile))
        for n in tre.traverse('postorder') :
            if not n.is_leaf() :
                if n.d not in dqs_sum :
                    dqs_sum[n.d] = [n.qs*num, num]
                else :
                    dqs_sum[n.d][0] += n.qs * num
                    dqs_sum[n.d][1] += num
    dqs_sum = { k: v[0] / v[1] for k, v in dqs_sum.items() }
    dqs_sum = { k: int(100*(v[0] - np.max(v[1:]))) for k, v in dqs_sum.items() }

    for n in tre.traverse('postorder') :
        if not n.is_leaf() :
            n.name = int(dqs_sum[n.d])
    tre.set_outgroup(tre.get_midpoint_outgroup())
    tre.unroot()
    with open(outfile+'.astral', 'wt') as fout :
        fout.write(tre.write(format=1)+'\n')
    return outfile+'.astral'

@click.command()
@click.option('-t', '--tree', help='supertree', required=True)
@click.option('-g', '--gene', help='genetrees', required=True)
@click.option('-d', '--dist', help='distance matrix', required=True)
@click.option('-o', '--outfile', help='output', required=True)
@click.option('-s', '--n_split', help='split gene trees into multiple parts for memory efficiency', default=1, type=int)
@click.option('-n', '--n_thread', help='threads to use', default=8, type=int)
def main(tree, gene, dist, outfile, n_split, n_thread) :
    iqr_tree = assign_iqr(outfile, tree, gene, n_split, n_thread)
    with uopen(dist) as fin, open(outfile+'.dist', 'wt') as fout :
        header = fin.readline().strip()
        fout.write('1\n\n%genome\n{0} 1000000\n'.format(header))
        for line in fin :
            fout.write(line)
    p = subprocess.Popen('{0} -i {2}.dist -o {2} -t {1}'.format(erable, iqr_tree, outfile).split(), stdout=subprocess.PIPE)
    p.communicate()
    os.unlink(iqr_tree)
    os.unlink(outfile+'.dist')
    tre = Tree('{0}.length.nwk'.format(outfile), format=1)
    for node in tre.iter_descendants() :
        if node.dist < 0 :
            sibling = [c for c in node.up.children if c != node]
            if len(sibling) :
                sibling[0].dist -= node.dist
            node.dist = 0

    tre.set_outgroup(tre.get_midpoint_outgroup())
    tre.unroot()
    with open(outfile, 'wt') as fout:
        fout.write(tre.write(format=1)+'\n')
    os.unlink(outfile+'.length.nwk')
    os.unlink(outfile+'.rates.txt')



if __name__ == '__main__' :
    main()