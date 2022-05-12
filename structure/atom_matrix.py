# coding=utf-8

import numpy as np
from gmx.chem.chem import Elements
from gmx.io.baseio import GMXDomains
from gmx.io.gro import readline_gro
from gmx.io.pdb import readline_pdb
from gmx.math.common import boxpara
from gmx.other.data_type import GMXDataType
from gmx.other.input_func import is_int
from gmx.structure.io_matrix import IOMatrix


class AtomMatrix(IOMatrix):
    def __init__(self):
        self.atomn = 0
        self.atoms = []
        self.elements = []
        self.molids = []
        self.molnms = []
        self.xyzs = np.array([], dtype=GMXDataType.REAL)
        self.velos = np.array([], dtype=GMXDataType.REAL)
        self.bonds = np.array([], dtype=GMXDataType.INT)

    def append(self, atomnm, element, molid, molnm, x, y, z, vx=0., vy=0., vz=0.):
        """Add atom record to this matrix"""
        self.atomn += 1
        self.atoms.append(atomnm)
        self.elements.append(element)
        self.molids.append(molid)
        self.molnms.append(molnm)
        self.xyzs = np.append(self.xyzs, [x, y, z])
        self.velos = np.append(self.velos, [vx, vy, vz])

    def clean(self):
        """Convert coordinate and velocity matrix to matrix with 3 columns"""
        self.xyzs = np.asarray(self.xyzs, dtype=GMXDataType.REAL).reshape(-1, 3)
        self.velos = np.asarray(self.velos, dtype=GMXDataType.REAL).reshape(-1, 3)
        self.bonds = np.asarray(self.bonds, dtype=GMXDataType.INT).reshape(-1, 2)

    def from_pdb(self, fpath):
        """Read atoms from PDB format file"""
        with open(fpath, 'r') as f:
            rown = 0
            molreindex = 0
            molkey = (0, '')
            for line in f.readlines():
                rown += 1
                if line[0:6] in ("ATOM  ", "HETATM"):
                    atominfo = readline_pdb(line, rown)
                    resn = atominfo[GMXDomains.RESID]
                    resnm = atominfo[GMXDomains.RESNM]
                    # Start a new molecule
                    if molkey != (resn, resnm):
                        molreindex += 1
                    self.append(
                        atominfo[GMXDomains.ATOMNM], atominfo[GMXDomains.ELEMENT],
                        molreindex, resnm,
                        *atominfo[GMXDomains.XYZ])
                    molkey = (resn, resnm)
            self.clean()

    def from_gro(self, fpath):
        """Read atoms from GRO format file"""
        with open(fpath, 'r') as f:
            f.readline()   # Skip title
            # Second line should be an integer equals to total atom numbers
            atomn = f.readline().strip()
            if is_int(atomn):
                atomn = int(atomn)
            else:
                raise ValueError('Unrecognized GRO format file. Expect ')
            rown = 2
            molreindex = 0
            molkey = (0, '')
            for line in f.readlines():
                rown += 1
                # Line 'atomn + 3' should be solvate box parameters
                if rown < atomn + 3:
                    atominfo = readline_gro(line, rown)
                    moln = atominfo[GMXDomains.MOLID]
                    molnm = atominfo[GMXDomains.MOLNM]
                    # Start a new molecule
                    if molkey != (moln, molnm):
                        molreindex += 1
                    atomnm = atominfo[GMXDomains.ATOMNM]
                    element = Elements.guess_element(atomnm)
                    self.append(
                        atomnm, element,
                        molreindex, molnm,
                        *atominfo[GMXDomains.XYZ],
                        *atominfo[GMXDomains.VXYZ])
                    molkey = (moln, molnm)
            self.clean()

    def to_pdb(self, fpath):
        """Write atoms to PDB format file"""
        outbuf = ""

        atomid = 0
        with open(fpath, 'w') as f:
            for atomnm, element, molnm, molid, xyz in zip(
                    self.atoms, self.elements, self.molnms, self.molids, self.xyzs):
                atomid += 1
                outbuf += IOMatrix.PDB_FORMAT.format(
                    atomid, atomnm, molnm, molid, *xyz, 1., 0., element)
                # End line with LF
                outbuf += '\n'
            
            if self.bonds.size != 0:
                connects = {}
                for b1, b2 in self.bonds:
                    if b1 in connects:
                        connects[b1].add(b2)
                    else:
                        connects[b1] = set([b2])
                    if b2 in connects:
                        connects[b2].add(b1)
                    else:
                        connects[b2] = set([b1])
                for ba, bc in connects.items():
                    bc = sorted(list(bc))
                    outbuf += f"CONECT{ba:>5d}" + ''.join([f'{a:>5d}' for a in bc]) + '\n'
            f.write(outbuf)
            f.write("END\n")

    def to_gro(self, fpath):
        """Write atoms to GRO format file"""
        outbuf = ""

        global_with_velo = True if np.linalg.norm(self.velos) > 1e-6 else False

        atomid = 0
        with open(fpath, 'w') as f:
            # Default title
            f.write("GROtesk MACabre and Sinister\n")
            # Atom numbers
            f.write("{0:>5}\n".format(self.atomn))
            for atomnm, molnm, molid, xyz, vxyz in zip(
                    self.atoms, self.molnms, self.molids, self.xyzs, self.velos):
                atomid += 1
                outbuf += IOMatrix.GRO_FORMAT.format(
                    molid, molnm, atomnm, atomid, *xyz)
                if global_with_velo:
                    outbuf += IOMatrix.GRO_VELOCITY.format(*vxyz)
                # End line with LF
                outbuf += '\n'
            f.write(outbuf)
            # Solvate box parameters
            f.write("{0:>10.5f}{1:>10.5f}{2:>10.5f}\n".format(
                *(boxpara(self.xyzs)[0])))
    
    def to_xyz(self, fpath):
        """Write atoms to XYZ format file"""
        outbuf = ""

        with open(fpath, 'w') as f:
            f.writelines(IOMatrix.XYZ_FORMAT.format(*xyz) + '\n' for xyz in self.xyzs)
            # End XYZ file with an extra empty line
            f.write('\n')
