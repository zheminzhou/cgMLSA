import os, sys, click, re, numpy as np
from ete3 import Tree, TreeNode
try :
    from .configure import externals, uopen
except :
    from configure import externals, uopen


def get_subgroup(tre, nodes, min_group, max_group) :
    tre.id, tre.tag = 0, '__G_0__'
    groups = [tre]
    min_groups = (min_group / (2 ** np.arange(np.log2(min_group).astype(int) - 1))).astype(int).tolist()
    if min_groups[-1] > 4 :
        min_groups.append(4)
    for min_group in min_groups :
        if np.max([len(g.get_leaves()) for g in groups]) <= max_group :
            break
        working = nodes
        nodes = []
        for dqs, node in working :
            root = node.get_tree_root()
            dNode = len(node.get_leaf_names())
            dTree = len(root.get_leaf_names()) + (root != tre)
            if dTree <= max_group :
                continue
            elif dNode + 1 < min_group or dTree - dNode + 1 < min_group :
                if dNode + 1 >= 4 and dTree - dNode + 1 >= 4 :
                    nodes.append([dqs, node])
                continue
            # split tree by node
            node.id = len(groups)
            node.tag = '__G_{0}__'.format(node.id)
            tagged_node = TreeNode(name=node.tag, dist=node.dist)
            parent = node.up
            parent.remove_child(node)
            node.up = None
            parent.add_child(tagged_node)
            tagged_node.up = parent
            groups.append(node)
    #groups = [ [ t.tag, t.dqs, t.get_leaf_names(), []] for t in groups ]
    roots = [None for t in groups]
    for t in groups :
        descs = [ int(re.findall(r'__G_(\d+)__', n)[0]) for n in t.get_leaf_names() if n.startswith('__G_') ]
        for leaf_id in descs :
            roots[leaf_id] = t.tag
    for t2, r in zip(groups, roots) :
        if r :
            n, p = TreeNode(dist=1., name=t2.name), TreeNode(name=r, dist=1.)
            children = t2.children[:]
            for c in children :
                t2.remove_child(c)
                n.add_child(c)
                c.up = n
            t2.add_child(n)
            t2.add_child(p)
            n.up, p.up = t2, t2
    return groups


@click.command()
@click.option('-i', '--tree', help='tree scored using ASTRAL', required=True)
@click.option('-o', '--output', help='sub-groups', required=True)
@click.option('-v', '--min_pp', help='min delta quartet support [default: 0.95]', default=0.95, type=float)
@click.option('-m', '--min_group', help='min size of sub-group [default: 100]', default=100, type=int)
@click.option('-x', '--max_group', help='max size of sub-group [default: 300]', default=300, type=int)
def main(tree, output, min_group, max_group, min_pp) :
    tre = Tree(tree, format=0, quoted_node_names=True)
    nodes = sorted([ [node.support, node] for node in tre.get_descendants('postorder') if node.support >= min_pp ], key=lambda n:-n[0])
    groups = get_subgroup(tre, nodes, min_group, max_group)

    with uopen(output, 'w') as fout :
        for t in groups :
            fout.write('{0}\t{1}\t{2:.2f}\t|\t{3}\n'.format(t.tag, len(t.get_leaves()), t.support, t.write(format=0)))

if __name__ == '__main__' :
    main()
