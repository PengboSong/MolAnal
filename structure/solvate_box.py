# coding=utf-8

import numpy as np
from gmx.math.common import boxpara
from gmx.other.data_type import GMXDataType


def unpack_args(*args, defv=0.):
    res = np.zeros(3, dtype=GMXDataType.REAL)
    if len(args) == 0:   # [defv, defv, defv]
        res.fill(defv)
    elif len(args) == 1:   # [x, x, x]
        x = args[0]
        res.fill(x)
    elif len(args) == 2:   # [x, x, z]
        x, z = args
        res[:] = [x, x, z]
    else:
        res[:] = args[:3]
    return res


class SolvateBox(object):
    def __init__(self):
        self.boxpara = np.zeros(3, dtype=GMXDataType.REAL)
        self.boxfactor = np.ones(3, dtype=GMXDataType.REAL)
        self.shift = np.zeros(3, dtype=GMXDataType.REAL)

    def set_boxpara(self, *args):
        self.boxpara[:] = unpack_args(*args)

    def set_boxfactor(self, *args):
        self.boxfactor[:] = unpack_args(*args, defv=1.)

    def clear_shift(self):
        self.shift.fill(0.)
    
    def need_shift(self):
        return np.linalg.norm(self.shift) > 1e-3

    def calc_boxpara(self, xyz):
        resizebox, shift = boxpara(xyz, extend=self.boxfactor)
        self.boxpara[:] = resizebox
        self.shift[:] = shift
        if self.need_shift():
            self.do_shift()
            self.clear_shift()

    def extend_box(self, x, y, z):
        boxvec = np.asarray([x, y, z], dtype=GMXDataType.REAL)
        self.boxpara += boxvec
        self.shift += .5 * boxvec
        if self.need_shift():
            self.do_shift()
            self.clear_shift()
    
    def extend_z(self, z):
        self.extend_box(x=0., y=0., z=z)

    def do_shift(self):
        """Shift system coordinate"""
        raise NotImplementedError('Shift system coordinate not supported.')
