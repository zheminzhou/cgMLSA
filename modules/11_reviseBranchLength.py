import sys, os, click, numpy as np
from ete3 import Tree
try :
    from .configure import uopen
except :
    from configure import uopen

@click.command()
@click.option('-t', '--tree', help='tree file', required=True)
@click.option('-d', '--dist', help='distance matrix', required=True)
@click.option('-o', '--outfile', help='output', required=True)
def main(tree, dist, outfile) :
    # load matrix
    with uopen(dist) as fin :
        headline = fin.readline()
        cnt = int(headline.strip())
        names = []
        distance = np.zeros([cnt, cnt], dtype=np.float32)
        for i, line in enumerate(fin) :
            part = line.strip().split()
            names.append(part[0])
            distance[i] = np.array(part[1:]).astype(np.float32)
        names = { n:i for i, n in enumerate(names) }
        dist_sum = np.sum(distance, 1)
    # load tree
    tre = Tree(tree, format=1)
    tre.set_outgroup(tre.get_midpoint_outgroup())
    for i, n in enumerate(tre.traverse()) :
        n.tag = i
        if n.is_leaf() :
            n.id = names[n.name]
    leaves = tre.get_leaves()
    while len(leaves) > 1 :
        #   get currently open merging
        merges = {}
        for leaf in leaves :
            leaf.id = leaf.id
            if leaf.up :
                if leaf.up.tag not in merges :
                    merges[leaf.up.tag] = [leaf]
                else :
                    merges[leaf.up.tag].append(leaf)
        merges = [ [pi, pj, None, None, None] for pairs in merges.values() \
                   if len(pairs) > 1 \
                   for i, pi in enumerate(pairs) for pj in pairs[:i] ]
        #   find minimum d'
        for m in merges :
            m[2] = (cnt-2)*distance[m[0].id, m[1].id] - dist_sum[m[0].id] - dist_sum[m[1].id]
        m = min(merges, key=lambda m:m[2])
        #   update branch lengths
        if cnt > 2 :
            m[3] = 0.5*distance[m[0].id, m[1].id] + 1./(2.*(cnt-2))*(dist_sum[m[0].id] - dist_sum[m[1].id])
            #if m[3] < 0 :
            #    m[3] = 0
            #elif m[3] > distance[m[0].id, m[1].id] :
            #    m[3] = distance[m[0].id, m[1].id]
            m[4] = distance[m[0].id, m[1].id] - m[3]
        else :
            m[3] = 0.5*distance[m[0].id, m[1].id]
            m[4] = 0.5*distance[m[0].id, m[1].id]

        #   update open nodes
        m[0].dist, m[1].dist = m[3], m[4]
        if min(m[3], m[4]) < -1 :
            print(m[0], m[0].dist, m[1], m[1].dist)
        new_node = m[0].up
        new_node.id = m[0].id
        distance[new_node.id, :] = 0.5*(distance[m[0].id, :] + distance[m[1].id, :] - distance[m[0].id, m[1].id])
        distance[:, new_node.id] = distance[new_node.id, :]
        distance[m[1].id, :] = 0
        distance[:, m[1].id] = 0
        dist_sum = np.sum(distance, 1)
        leaves = [ leaf for leaf in leaves if leaf not in m[:2] ] + [new_node]
        cnt = len(leaves)
    with uopen(outfile, 'wt') as fout :
        fout.write(tre.write(format=1)+'\n')

if __name__ == '__main__' :
    main()