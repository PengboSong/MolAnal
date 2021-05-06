# coding=utf-8

import numpy as np
from gmx.chem.chem import Elements
from gmx.io.gro import readline_gro
from gmx.io.pdb import readline_pdb
from gmx.math.common import boxpara
from gmx.other.data_type import GMXDataType
from gmx.other.input_func import is_int
from gmx.other.mol_order import MolOrder

from io_matrix import IOMatrix


class SingleMol(object):
    def __init__(self, molid, molnm):
        self.i = molid
        self.nm = molnm
        self.start_atomid = 0
        self.atomn = 0
        self.atoms = []
        self.elements = []
        self.xyzs = np.array([], dtype=GMXDataType.REAL)
        self.velos = np.array([], dtype=GMXDataType.REAL)

    def append(self, atomnm, element, x, y, z, vx=0., vy=0., vz=0.):
        """Add atom record to this molecule"""
        self.atomn += 1
        self.atoms.append(atomnm)
        self.elements.append(element)
        self.xyzs = np.append(self.xyzs, [x, y, z])
        self.velos = np.append(self.velos, [vx, vy, vz])

    def clean(self):
        """Convert coordinate and velocity matrix to matrix with 3 columns"""
        self.xyzs = np.asarray(self.xyzs, dtype=GMXDataType.REAL).reshape(-1, 3)
        self.velos = np.asarray(self.velos, dtype=GMXDataType.REAL).reshape(-1, 3)

    def to_pdb(self):
        """Write molecule to PDB format lines"""
        outbuf = ""

        atomid = self.start_atomid
        for atomnm, element, xyz in zip(self.atoms, self.elements, self.xyzs):
            atomid += 1
            outbuf += IOMatrix.PDB_FORMAT.format(
                atomid, atomnm, self.nm, self.i, *xyz, 1., 0., element)
        return outbuf

    def to_gro(self):
        """Write molecule to GRO format lines"""
        outbuf = ""

        with_velo = True if np.linalg.norm(self.velos) > 1e-6 else False

        atomid = self.start_atomid
        for atomnm, xyz, vxyz in zip(self.atoms, self.xyzs, self.velos):
            atomid += 1
            outbuf += IOMatrix.GRO_FORMAT.format(
                self.i, self.nm, atomnm, atomid, *xyz)
            if with_velo:
                outbuf += IOMatrix.GRO_VELOCITY.format(*vxyz)
        return outbuf


class MolMatrix(IOMatrix):
    def __init__(self):
        self.mols = []

    def from_pdb(self, fpath):
        """Read molecules from PDB format file"""
        with open(fpath, 'r') as f:
            rown = 0
            molreindex = 0
            molkey = (0, '')
            molmat = None
            for line in f.readlines():
                rown += 1
                if line[0:6] in ("ATOM  ", "HETATM"):
                    atominfo = readline_pdb(line, rown)
                    resn = atominfo["residue_sequence_number"]
                    resnm = atominfo["residual_name"]
                    # Start a new molecule
                    if molkey != (resn, resnm):
                        if molmat:
                            molmat.clean()
                            self.mols.append(molmat)
                        molreindex += 1
                        molmat = SingleMol(molreindex, resnm)
                    molmat.append(
                        atominfo["atom_name"], atominfo["element_symbol"],
                        *atominfo["coordinate"])
                    molkey = (resn, resnm)
            molmat.clean()
            self.mols.append(molmat)
            self.resort(MolOrder.KEEP_ORDER)

    def from_gro(self, fpath):
        """Read molecules from GRO format file"""
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
            molmat = None
            for line in f.readlines():
                rown += 1
                # Line 'atomn + 3' should be solvate box parameters
                if rown < atomn + 3:
                    atominfo = readline_gro(line, rown)
                    moln = atominfo["mol_id"]
                    molnm = atominfo["mol_name"]
                    # Start a new molecule
                    if molkey != (moln, molnm):
                        if molmat:
                            molmat.clean()
                            self.mols.append(molmat)
                        molreindex += 1
                        molmat = SingleMol(molreindex, molnm)
                    atomnm = atominfo["atom_name"]
                    element = Elements.guess_element(atomnm)
                    molmat.append(
                        atomnm, element,
                        *atominfo["coordinate"],
                        *atominfo["velocity"])
                    molkey = (moln, molnm)
            molmat.clean()
            self.mols.append(molmat)
            self.resort(MolOrder.KEEP_ORDER)

    def to_pdb(self, fpath):
        """Write molecules to PDB format file"""
        with open(fpath, 'w') as f:
            f.writelines(molmat.to_pdb() + '\n' for molmat in self.mols)
            f.write("END\n")

    def to_gro(self, fpath):
        """Write molecules to GRO format file"""
        system_xyzs = np.vstack([molmat.xyzs for molmat in self.mols])
        total_atomn = system_xyzs.shape[0]
        with open(fpath, 'w') as f:
            # Default title
            f.write("GROtesk MACabre and Sinister\n")
            # Atom numbers
            f.write("{0:>5}\n".format(total_atomn))
            f.writelines(molmat.to_gro() + '\n' for molmat in self.mols)
            # Solvate box parameters
            f.write("{0:>10.5f}{1:>10.5f}{2:>10.5f}\n".format(
                *boxpara(system_xyzs)))

    def resort(self, ordertype):
        if ordertype not in MolOrder:
            ordertype = MolOrder.KEEP_ORDER

        count_atomn = 0
        if ordertype == MolOrder.KEEP_ORDER:
            for mol in self.mols:
                mol.start_atomid = count_atomn
                count_atomn += mol.atomn
        elif ordertype in [MolOrder.MOL_ORDER, MolOrder.MOL_ORDER_ALPHA]:
            molreindex = 0
            mol_groups = self.clustering()
            mol_types = mol_groups.keys() if ordertype == MolOrder.MOL_ORDER \
                                          else sorted(mol_groups.keys())
            for moltype in mol_types:
                for j in mol_gsroups[moltype]:
                    thismol = self.mols[j - 1]
                    assert thismol.i == j
                    molreindex += 1
                    thismol.i = molreindex
                    thismol.start_atomid = count_atomn
                    count_atomn += mol.atomn
            self.mols.sort(key=lambda mol: mol.i)

    def clustering(self):
        """Clustering molecules into groups based on molecular types"""
        groups = {}
        for mol in self.mols:
            molnm = mol.nm
            if molnm in groups:
                groups[molnm].append(mol.i)
            else:
                groups[molnm] = [mol.i]
        # Sort molecules by ID in ascending order
        for molnm in groups:
            groups[molnm].sort()
        return groups
