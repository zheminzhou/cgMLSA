import os, sys, numpy as np, pandas as pd, click, time
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, OPTICS, cluster_optics_dbscan, DBSCAN
from sklearn.metrics import silhouette_score, silhouette_samples


def PCA_Kmeans(prefix, names, dist) :
    if len(names) <= 1 :
        return [[prefix, names, dist]]
    model = PCA(n_components=min(len(names)-1, 10), svd_solver='randomized', random_state=43, copy=True)
    components = model.fit_transform(dist)

    main_idx = model.explained_variance_ratio_/model.explained_variance_ratio_[0] >= 0.01
    c = components[:, main_idx]
    res = [0.5, 1, np.zeros(c.shape[0], dtype=int)]

    plt.figure(figsize=(8,8))
    for n in np.exp((np.log(10.) - np.log(0.01))*range(20)/20.+np.log(0.01)) :
        #kmeans = KMeans(n_clusters=n, n_init=100, init='k-means++').fit(c)
        clust = DBSCAN(eps=n).fit(c)
        labels = clust.labels_
        if np.max(labels) > 1 :
            ss = silhouette_score(dist[labels>=0][:, labels>=0], labels=labels[labels>=0], metric='precomputed')
        else :
            ss = 0.
        print('{0} {1} {2}'.format(prefix, n, ss))

        plt.subplot(311)
        if components.shape[1] > 1 :
            plt.scatter(components[labels<0, 0], components[labels<0, 1], c='grey')
            plt.scatter(components[labels>=0, 0], components[labels>=0, 1], c=labels[labels>=0])
        if components.shape[1] > 2 :
            plt.subplot(312)
            plt.scatter(components[labels<0, 0], components[labels<0, 2], c='grey')
            plt.scatter(components[labels>=0, 0], components[labels>=0, 2], c=labels[labels>=0])
            plt.subplot(313)
            plt.scatter(components[labels<0, 1], components[labels<0, 2], c='grey')
            plt.scatter(components[labels>=0, 1], components[labels>=0, 2], c=labels[labels>=0])
        plt.show()
        if ss > res[0] :
            res = [ss, n, labels]

    outputs = []
    for i in range( np.max(res[2])+1 ) :
        outputs.append([ '{0}.{1}'.format(prefix, i+1), names[res[2] == i], dist[res[2] == i][:, res[2] == i] ])
    return outputs


@click.command()
@click.option('-d', '--distfile', help='dist file', default='')
@click.option('-z', '--npzfile', help='npz file', default='')
def main(distfile, npzfile) :
    if npzfile :
        zdata = np.load(npzfile)
        names, dist = zdata['names'], zdata['dist']
    else :
        with open(distfile) as fin :
            fin.readline()
            dist = pd.read_csv(fin, sep=' ', dtype=str, header=None).values
            names, dist = dist.T[0].astype(str), dist[:, 1:].astype(np.float32)
    groups = [['0', names, dist]]
    new_groups = []
    for ite in range(5) :
        for group in groups :
            ng = PCA_Kmeans(*group)
            new_groups.extend(ng)
        groups = new_groups
    return new_groups


if __name__ == '__main__' :
    main()
