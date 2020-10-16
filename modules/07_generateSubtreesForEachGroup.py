import os, sys, click, re, numpy as np
from ete3 import Tree, TreeNode
from copy import deepcopy
from multiprocessing import Pool
try :
    from .configure import externals, uopen
except :
    from configure import externals, uopen

def pick(subset, g1, g2, num) :
    x = list(subset[g1] - {g2})
    chosen = []
    while num > len(x) :
        chosen.extend(np.random.permutation(x))
        num -= len(x)
    chosen.extend(np.random.choice(x, num, replace=False))

    secondary_choice = {}
    for c in chosen :
        if c in subset :
            secondary_choice[c] = secondary_choice.get(c, 0)+1
    if len(secondary_choice) :
        for g, cnt in secondary_choice.items() :
            secondary_choice[g] = pick(subset, g, g1, cnt)
        for i, c in enumerate(chosen) :
            if c in subset :
                chosen[i] = secondary_choice[c][-1]
                secondary_choice[c] = secondary_choice[c][:-1]
    return chosen

def get_a_subtree(tre, tips) :
    nodes = []
    n_tip = 0
    for node in tre.traverse('postorder'):
        if node.is_leaf():
            if node.name in tips:
                node.d = len(nodes)
                nodes.append(TreeNode(name=tips[node.name], dist=node.dist))
                n_tip += 1
            else:
                node.d = -1
        else:
            desc = [c.d for c in node.children if c.d >= 0]
            if len(desc) == 0:
                node.d = -1
            elif len(desc) == 1:
                node.d = desc[0]
                nodes[desc[0]].dist += node.dist
            else:
                node.d = len(nodes)
                nodes.append(TreeNode(name=node.name, dist=node.dist))
                for id in desc:
                    nodes[-1].add_child(nodes[id])
                    nodes[id].up = nodes[-1]
    if n_tip >= 4 :
        return nodes[-1]
    else :
        return None

def get_subtrees(data) :
    prefix, tree, subsets, repeats = data
    tre = Tree(tree, format=1)
    leaves = {n: 1 for n in tre.get_leaf_names()}

    selected = { grp:{ tip:tip for tip in tips if tip in subsets or tip in leaves } for grp, tips in subsets.items() }
    for grp, tips in selected.items() :
        for tip in tips :
            if tip in subsets :
                tips[tip] = [ t if t in leaves else None for t in pick(subsets, tip, grp, repeats)]

    results = {}
    for grp, mtips in selected.items() :
        results[grp] = []
        if len(mtips) < 4 :
            continue
        possible_tips = {}
        for tip, tgrp in mtips.items() :
            if tip in subsets :
                possible_tips.update( { t:t for t in tgrp } )
            else :
                possible_tips[tip] = tip
        subtree1 = get_a_subtree(tre, possible_tips)
        if subtree1 :
            for i in range(repeats) :
                tips = {}
                for tip, tgrp in mtips.items() :
                    if tip in subsets :
                        tips[tgrp[i]] = tip
                    else :
                        tips[tip] = tip
                subtree = get_a_subtree(subtree1, tips)
                if subtree :
                    subtree.unroot()
                    results[grp].append(subtree.write(format=1))
    return results

@click.command()
@click.option('-i', '--trees', help='gene trees', required=True)
@click.option('-g', '--groups', help='groups', required=True)
@click.option('-p', '--prefix', help='prefix', required=True)
@click.option('-r', '--repeats', help='number of sampling per tree', default=10)
def main(trees, groups, prefix, repeats) :
    subsets = {}
    with uopen(groups) as fin :
        for line in fin :
            part = line.strip().split('\t')
            subsets[part[0]] = set(Tree(part[4], format=1).get_leaf_names()) #set(part[4:])

    with uopen(trees) as fin :
        trees = fin.readlines()

    results = { grp:[] for grp in subsets }
    tree_id = 0
    for ts in pool.imap_unordered(get_subtrees, [ [prefix, tree, subsets, repeats] for tree in trees ]) :
        tree_id += 1
        if tree_id % 100 == 0 :
            sys.stderr.write('Sub-grouped {0} gene trees.\n'.format(tree_id))
        for grp, t in ts.items() :
            results[grp].extend(t)
    for grp, ts in results.items() :
        with open('{0}.subset.{1}'.format(prefix, grp), 'w') as fout:
            fout.write('\n'.join(ts)+'\n')

pool = Pool(10)
if __name__ == '__main__' :
    main()