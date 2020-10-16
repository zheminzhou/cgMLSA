import click, numpy as np, sys

def readMat(dist) :
    names, distances = [], []
    with open(dist) as fin :
        fin.readline()
        for line in fin :
            part = line.strip().split()
            names.append(part[0])
            distances.append([ float(d) for d in part[1:] ])
    return names, distances

@click.command()
@click.argument('dist1')
@click.argument('dist2')
def main(dist1, dist2) :
    n1, d1 = readMat(dist1)
    n2, d2 = readMat(dist2)
    idx = { n:i for i, n in enumerate(n1) }
    idx1 = [ idx[n] for n in n2 if n in idx ]
    idx2 = [ i for i, n in enumerate(n2) if n in idx ]
    n1 = np.array(n1)[np.array(idx1)]
    n2 = np.array(n2)[np.array(idx2)]
    d1 = np.array(d1)[idx1, :][:, idx1]
    d2 = np.array(d2)[idx2, :][:, idx2]
    from _collections import defaultdict
    for i, nn1 in enumerate(n1) :
        for j, nn2 in enumerate(n2[:i]) :
            print(nn1, nn2, d1[i, j], d2[i, j])

if __name__ == '__main__' :
    main()