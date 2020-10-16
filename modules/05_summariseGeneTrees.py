
import sys, os, click, subprocess, numpy as np, re
from ete3 import Tree
try :
    from .configure import externals
except :
    from configure import externals

astrid = externals['astrid']
astral = externals['astral']

def assign_qs(tree) :
    tre = Tree(tree, format=1, quoted_node_names=True)
    tre.qs = np.array([1., 0., 0.])
    for node in tre.iter_descendants('postorder') :
        if not node.is_leaf() :
            qs = re.findall('\[q1=(.+);q2=(.+);q3=(.+)\]', node.name)
            try :
                node.qs = np.array([ float(q) for q in qs[0]])
            except Exception :
                node.qs = np.zeros(3, dtype=float)
    for node in tre.traverse('postorder') :
        if node.is_leaf() :
            node.d = (node.name,)
        else :
            node.d = tuple(sorted([ n for c in node.children for n in c.d ]))
    return tre

@click.command()
@click.option('-i', '--input', help='input file')
@click.option('-o', '--output', help='output file')
@click.option('-s', '--split', help='split gene trees into multi-parts in ASTRAL step for memory efficiency', default=1, type=int)
@click.option('-t', '--n_thread', help='num of threads in use', default=10, type=int)

def main(input, output, split, n_thread) :
    # run astrid
    p = subprocess.Popen('{0} {3} -i {1} -o {2}.astrid -c {2}.atid'.format(astrid, input, output).split(), stdout=subprocess.PIPE)
    p.communicate()
    if p.returncode > 0 :
        sys.exit(1)

    if split <= 1 :
        geneTrees = [[1, input]]
    else :
        geneTrees = []
        with open(input) as fin :
            trees = fin.readlines()
        for ite in np.arange(split) :
            s_file = output+'.split.{0}'.format(ite)
            with open(s_file, 'wt') as fout :
                subtrees = trees[ite::split]
                fout.write('\n'.join(subtrees))
            geneTrees.append([len(subtrees), s_file])

    # run astral
    dqs_sum = {}
    tre = None
    for i, (num, tfile) in enumerate(geneTrees) :
        p = subprocess.Popen('{0} -q {1}.astrid -i {2} -t 8 -o {1}.astral.{3} -T {4}'.format(astral, output, tfile, i, n_thread).split(), \
                             stdout=subprocess.PIPE)
        p.communicate()
        if p.returncode > 0 :
            sys.exit(1)
        tre = assign_qs('{0}.astral.{1}'.format(output, i))
        os.unlink(tfile)
        os.unlink('{0}.astral.{1}'.format(output, i))
        for n in tre.traverse('postorder') :
            if not n.is_leaf() :
                if n.d not in dqs_sum :
                    dqs_sum[n.d] = [n.qs*num, num]
                else :
                    dqs_sum[n.d][0] += n.qs * num
                    dqs_sum[n.d][1] += num
    dqs_sum = { k: v[0] / v[1] for k, v in dqs_sum.items() }
    dqs_sum = { k: v[0] - np.max(v[1:]) for k, v in dqs_sum.items() }
    for n in tre.traverse('postorder') :
        if not n.is_leaf() :
            n.name = '{:.3f}'.format(dqs_sum[n.d])
    with open(output, 'wt') as fout :
        fout.write(tre.write(format=1)+'\n')

if __name__ == '__main__' :
  main()