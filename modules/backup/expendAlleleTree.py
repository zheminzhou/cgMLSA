import click, pandas as pd, re, numpy as np
from ete3 import Tree
from configure import readFasta, uopen, rc


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


def delTip(tre, tip):
    parent = tip.up
    if parent is None:
        return delNode(tre, tip)
    parent.remove_child(tip)
    tip.up = None
    if len(parent.children) < 2:
        tre = delNode(tre, parent)
    return tre


def getNumber(s, default=0):
    try:
        return float(s)
    except:
        return default


@click.command()
@click.option('-p', '--profile', help='profile', required=True)
@click.option('-f', '--fasta', help='fasta file', default=None)
@click.argument('treefile', nargs=1)
def main(profile, fasta, treefile):
    data = pd.read_csv(profile, delimiter='\t', header=0, dtype=str)
    loci = data.columns
    data = data.values
    genomes = data.T[0]
    tre = Tree(treefile, format=1)
    n = tre.get_leaf_names()[0].strip("'").rsplit('_', 1)[0]
    locId = np.where(loci == n)[0][0]

    tags = {}
    for g, v in zip(genomes, data[:, locId]):
        tag = '{0}_{1}'.format(n, v)
        if tag not in tags:
            tags[tag] = [g]
        else:
            tags[tag].append(g)
    labels = {}
    if fasta:
        seqs = readFasta(fasta)
        seq_len = len(list(seqs.values())[0])
        min_dist = 1. / 4. / seq_len
    else:
        min_dist = 1e-6

    max_node_support = 0.
    for node in tre.iter_descendants():
        if not node.is_leaf():
            s = getNumber(node.name)
            if s > max_node_support:
                max_node_support = s

    max_node_support = 100 if max_node_support > 1 else 1
    for node in tre.get_descendants('postorder'):
        if node.is_leaf():
            n = node.name.strip("'").split()[0]
            if n in tags:
                node.name = '{{{0}}}'.format(n)
                genomes = tags[n]
                if len(genomes) > 1:
                    if node.dist < min_dist:
                        labels[n] = '{0}'.format(','.join(['{0}:{1}'.format(g, node.dist) for g in genomes]))
                        node.dist = -1000
                    else:
                        labels[n] = '({0}){1}'.format(','.join(['{0}:0'.format(g) for g in genomes]), max_node_support)
                    node.name = '{{{0}}}'.format(n)
                else:
                    node.name = genomes[0]
            else:
                tre = delTip(tre, node)
        elif node.dist < min_dist:  # or getNumber(node.name, max_node_support) < 0.5*max_node_support :
            tre = delNode(tre, node)

    print(re.sub(r':-\d+\.*\d+', '', tre.write(format=1).format(**labels)))


if __name__ == '__main__':
    main()