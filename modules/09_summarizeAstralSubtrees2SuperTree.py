import os, glob, click, re
from ete3 import Tree

@click.command()
@click.option('-o', 'outfile', help='generated super-tree')
@click.option('-d', '--folder', help='folder storing gene trees', default='.')
@click.option('-i', '--input_pattern', help='Default: *.__G_*__.astral', default='*.__G_*__.astral')
def main(folder, input_pattern, outfile) :
    subtrees = glob.glob(os.path.join(folder, input_pattern))
    trees = [None] * len(subtrees)
    for treefile in subtrees :
        grp = re.findall('__G_(\d+)__', treefile)[0]
        tre = Tree(treefile, format=1, quoted_node_names=True)
        trees[int(grp)] = tre

    combined = {'__G_0__': 1}
    while True :
        descs = [ n for n in trees[0].get_leaves() if n.name.startswith('__G_') ]
        if len(descs) == 0 :
            break
        for linker1 in descs :
            combined[linker1.name] = 1
            id = int(re.findall('__G_(\d+)__', linker1.name)[0])
            tre2 = trees[id]
            linker2 = [ n for n in tre2.get_leaves() if n.name in combined ][0]
            tre2.set_outgroup(linker2)
            tre2 = [ c for c in tre2.children if c != linker2 ][0]
            parent = linker1.up
            parent.remove_child(linker1)
            parent.add_child(tre2)
            tre2.up = parent
            tre2.dist = 1.

    for n in trees[0].traverse() :
        if not n.is_leaf() :
            n.name = ''
    with open(outfile, 'wt') as fout :
        fout.write(trees[0].write(format=1) + '\n')

if __name__ == '__main__' :
    main()
