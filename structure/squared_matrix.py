# coding=utf-8

import numpy as np
from gmx.other.data_type import GMXDataType


class SquaredMatrix(object):
    def __init__(self, n, dtype=GMXDataType.REAL):
        self._N = n
        self._sqmat = np.zeros((n, n), dtype=dtype)
    
    def __getitem__(self, i):
        return self._sqmat[i]
    
    def __setitem__(self, i, val):
        self._sqmat[i] = val

    def copy(self):
        return self._sqmat.copy()

    def upper_half(self):
        """Return upper half of the squared matrix"""
        upper = np.array([])
        for i in range(self._N):
            upper = np.append(upper, self._sqmat[i, i+1:])
        return upper

    def elements(self):
        """Put every matrix elements into a seq with its i(row id) and j(column id)"""
        pairseq = []
        for i in range(self._N):
            for j in range(i+1, self._N):
                pairseq.append((i, j, self._sqmat[i, j]))
        return pairseq
