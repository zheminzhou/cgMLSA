
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
                if np.any(np.isnan(node.qs)) :
                    node.qs = np.zeros(3, dtype=float)
            except Exception :
                node.qs = np.zeros(3, dtype=float)
    for node in tre.traverse('postorder') :
        if node.is_leaf() :
            node.d = (node.name,)
        else :
            node.d = tuple(sorted([ n for c in node.children for n in c.d ]))
    return tre

@click.command()
@click.argument('trees', nargs=-1)
def main(trees) :
    dqs_sum = {}
    tre = None
    for i, tfile in enumerate(trees) :
        print(i)
        tre = assign_qs(tfile)
        for n in tre.traverse('postorder') :
            if not n.is_leaf() :
                if n.d not in dqs_sum :
                    dqs_sum[n.d] = [n.qs, 1]
                else :
                    dqs_sum[n.d][0] += n.qs
                    dqs_sum[n.d][1] += 1
    dqs_sum = { k: v[0] / v[1] for k, v in dqs_sum.items() }
    dqs_sum = { k: v[0] - np.max(v[1:]) for k, v in dqs_sum.items() }
    for n in tre.traverse('postorder') :
        if not n.is_leaf() :
            n.name = '{:.3f}'.format(dqs_sum[n.d])
    print(tre.write(format=1)+'\n')

if __name__ == '__main__' :
  main()