import click, numpy as np
from ete3 import Tree
from _collections import defaultdict

def delNode(tre, node):
    dist = node.dist
    parent = node.up
    if parent is None:
        if tre != node:
            return tre
        parent = node.children[0]
        dist = parent.dist
        node.remove_child(parent)
        parent.up = None
        tre = parent
    else:
        parent.remove_child(node)

    for child in node.children:
        parent.add_child(child)
        child.up = parent
        child.dist += dist
    return tre

def getBipartition(tre, names) :
    bipartition = {}
    for node in tre.traverse('postorder') :
        if node.is_leaf() :
            node.descendents = np.zeros(len(names), dtype=int)
            if node.name in names :
                node.descendents[names[node.name]] = 1
        else :
            node.descendents = np.sum([c.descendents for c in node.children], 0)
        bp = node.descendents if node.descendents[0] == 0 else 1 - node.descendents
        bipartition[tuple(bp.tolist())] = 1
    return bipartition

@click.command()
@click.argument('trees', nargs=-1)
def main(trees) :
    names = {}
    bipartition = defaultdict(int)
    for treefile in reversed(trees) :
        tre = Tree(treefile, format=1)
        if len(names) == 0 :
            names = {n:i for i, n in enumerate(sorted(tre.get_leaf_names()))}
        bipart = getBipartition(tre, names)
        for bp in bipart :
            bipartition[bp] += 1
    bipartition = { k:v for k,v in bipartition.items() if v >= len(trees) }
    for node in tre.traverse('postorder') :
        bp = node.descendents if node.descendents[0] == 0 else 1 - node.descendents
        if tuple(bp.tolist()) not in bipartition :
            tre = delNode(tre, node)
    print(tre.write(format=1))
if __name__ == '__main__' :
    main()