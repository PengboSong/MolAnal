# coding=utf-8

import numpy as np
from gmx.io.gro import readline_gro
from gmx.io.pdb import readline_pdb
from gmx.other.data_type import GMXDataType
from gmx.other.input_func import is_int

from io_matrix import IOMatrix


class CoordMatrix(IOMatrix):
    def __init__(self):
        self.atomn = 0
        self.xyzs = np.array([], dtype=GMXDataType.REAL)

    def append(self, x, y, z):
        """Add atom coordinate to this matrix"""
        self.atomn += 1
        self.xyzs = np.append(self.xyzs, [x, y, z])

    def clean(self):
        """Convert coordinate matrix to matrix with 3 columns"""
        self.xyzs = np.asarray(self.xyzs, dtype=GMXDataType.REAL).reshape(-1, 3)
    
    def from_matrix(self, xyz):
        """Pack a row matrix with 3 columns as coordinate matrix"""
        if xyz.ndim == 2 and xyz.shape[1] == 3:
            self.atomn = xyz.shape[0]
            self.xyzs = xyz
        else:
            print("[WARNING] Input matrix can not be packed into a coordinate matrix.")

    def from_pdb(self, fpath):
        """Read coordinate from PDB format file"""
        with open(fpath, 'r') as f:
            rown = 0
            for line in f.readlines():
                rown += 1
                if line[0:6] in ("ATOM  ", "HETATM"):
                    xyz = readline_pdb(line, rown)["coordinate"]
                    self.append(*xyz)
            self.clean()

    def from_gro(self, fpath):
        """Read coordinate from GRO format file"""
        with open(fpath, 'r') as f:
            f.readline()   # Skip title
            # Second line should be an integer equals to total atom numbers
            atomn = f.readline().strip()
            if is_int(atomn):
                atomn = int(atomn)
            else:
                raise ValueError('Unrecognized GRO format file. Expect ')
            rown = 2
            for line in f.readlines():
                rown += 1
                # Line 'atomn + 3' should be solvate box parameters
                if rown < atomn + 3:
                    xyz = readline_gro(line, rown)["coordinate"]
                    self.append(*xyz)
            self.clean()

    def to_xyz(self, fpath):
        """Write coordinate to XYZ format file"""
        outbuf = ""

        with open(fpath, 'w') as f:
            for xyz in self.xyzs:
                outbuf += IOMatrix.XYZ_FORMAT.format(*xyz)
                # End line with LF
                outbuf += '\n'
            # End XYZ file with an extra empty line
            f.write(outbuf + '\n')
