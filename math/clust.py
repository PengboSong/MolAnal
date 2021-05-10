# coding=utf-8

from gmx.other.data_type import GMXDataType
from gmx.structure.squared_matrix import SquaredMatrix


class Cluster(object):
    """Clustering N particles using given distance matrix.

    Attrs:
        clusterid: For each of N particles, list which cluster it is assigned.
        cluster: List all clusters and their members.
    """

    def __init__(self, n):
        self._N = n
        self._distmat = SquaredMatrix(n=n, dtype=GMXDataType.REAL)
        self.clusterid = list(range(n))
        self.cluster = {}

    def set_distance(self, sqmat):
        """Set distance matrix from outside source"""
        assert self._N == sqmat._N, "Source matrix should have the same size as this one."
        self._distmat.copy_from(sqmat)

    def link(self, cutoff):
        """Link clusters using given cutoff.
        
        If two particles have a distance less than cutoff, then these two
        particles are clustered into one group. Generated clustered group
        ids may be discontinuous, and clustered groups overview is not
        generated.
        """
        changed_status = True
        while (changed_status):
            changed_status = False
            for x in range(1, self._N):
                for y in range(x + 1, self._N):
                    if self._distmat[x, y] > cutoff:
                        continue
                    else:
                        xclust = self.clusterid[x]
                        yclust = self.clusterid[y]
                        diff = xclust - yclust
                        if diff != 0:
                            changed_status = True
                            if diff > 0:
                                self.clusterid[y].update(xclust)
                            else:
                                self.clusterid[x].update(yclust)

    def renum(self):
        """Renumber clustered group ids and generate groups overview"""
        old_clusterid = list(set(self.clusterid))
        new_clusterid = list(range(1, len(old_clusterid) + 1))
        cluster_reindex = dict(zip(old_clusterid, new_clusterid))
        for i, cid in enumerate(self.clusterid):
            ncid = cluster_reindex[cid]
            self.cluster[i] = ncid
            if ncid in self.cluster:
                self.cluster[ncid].append(i)
            else:
                self.cluster[ncid] = [i]
