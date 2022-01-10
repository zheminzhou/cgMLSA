
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
    for node in tre.traverse('postorder') :
        if node.is_leaf() :
            node.d = (node.name,)
        else :
            node.d = tuple(sorted([ n for c in node.children for n in c.d ]))
    return tre

@click.command()
@click.option('-i', '--input', help='input file')
@click.option('-o', '--output', help='output file')
@click.option('-t', '--n_thread', help='num of threads in use', default=10, type=int)
def main(input, output, n_thread) :
    # run astrid
    p = subprocess.Popen('{0} -i {1} -o {2}.astrid -c {2}.atid'.format(astrid, input, output).split(), stdout=subprocess.PIPE)
    p.communicate()
    if p.returncode > 0 :
        sys.exit(1)

    p = subprocess.Popen('{0} -q {1}.astrid -i {2} -t 3 -o {1} -T {3}'.format(astral, output, input, n_thread).split(), \
                         stdout=subprocess.PIPE)
    p.communicate()

if __name__ == '__main__' :
  main()
