# coding=utf-8

import numpy as np

from gmx.io.filetype import valid_filetype, valid_onlytype
from gmx.io.pdb import readline_pdb
from gmx.io.gro import readline_gro


class MolMatrix(object):
    def __init__(self, molid, molnm):
        self.i = molid
        self.nm = molnm
        self.atoms = []
        self.elements = []
        self.xyzs = []
        self.velos = []

        # PDB template line
        self.pdbformat = ("HETATM{0:>5} {1:<4} {2:>3}  {3:>4}    {4:>8.3f}{5:>8.3f}{6:>8.3f}"
                          "{7:>6.2f}{8:>6.2f}          {9:>3}\n")
        # GRO template line
        self.groformat = ""

    def append(self, atomnm, x, y, z, vx=0., vy=0., vz=0.):
        """Add atom record to this molecule."""
        self.atoms.append(atomnm)
        self.elements.append(atomnm.upper())
        self.xyzs.extend([x, y, z])
        self.velos.extend([vx, vy, vz])

    def clean(self):
        """Convert coordinate and velocity matrix to matrix with 3 columns."""
        self.xyzs = np.asarray(self.xyzs, dtype="float64").reshape(-1, 3)
        self.velos = np.asarray(self.velos, dtype="float64").reshape(-1, 3)

    def to_pdb(self):
        """Write molecule to PDB format lines"""
        outbuf = ""

        atomid = 0
        for atomnm, element, xyz in zip(self.atoms, self.elements, self.xyzs):
            atomid += 1
            outbuf += self.pdbformat.format(
                atomid, atomnm, self.nm, self.i, *xyz, 1., 0., element)
        return outbuf
    
    def to_gro(self):
        """Write molecule to GRO format lines"""
        pass #TODO

class MolModel(object):
    def __init__(self):
        self.mols = []

    def from_pdb(self, fpath):
        """Read molecules from PDB format file."""
        with open(fpath, 'r') as f:
            rown = 0
            molreindex = 0
            molkey = (0, '')
            molmat = None
            for line in f.readlines():
                if line[0:6] in ("ATOM  ", "HETATM"):
                    atominfo = readline_pdb(line, rown)
                    resn = atominfo["residue_sequence_number"]
                    resnm = atominfo["residual_name"]
                    if molkey != (resn, resnm):
                        if molmat:
                            molmat.clean()
                            self.mols.append(molmat)
                        molreindex += 1
                        molmat = MolMatrix(molreindex, resnm)
                    molmat.append(
                        atominfo["atom_name"], *atominfo["coordinate"])
                    molkey = (resn, resnm)
            molmat.clean()
            self.mols.append(molmat)
    
    def from_gro(self, fpath):
        """Read molecules from GRO format file."""
        pass #TODO
    
    def to_pdb(self, fpath):
        """Write molecules to PDB format file."""
        with open(fpath, 'w') as f:
            f.writelines(molmat.to_pdb() for molmat in self.mols)
            f.write("END\n")
    
    def to_gro(self, fpath):
        """Write molecules to GRO format file."""
        pass #TODO

    def from_file(self, fpath):
        """Read molecules from file."""
        if valid_filetype(fpath, filetype="pdb"):
            self.from_pdb(fpath)
        elif valid_filetype(fpath, filetype="gro"):
            self.from_gro(fpath)
    
    def to_file(self, fpath):
        """Write molecules to file."""
        if valid_onlytype(fpath, filetype="pdb"):
            self.to_pdb(fpath)
        elif valid_onlytype(fpath, filetype="gro"):
            self.to_gro(fpath)
    
    def clustering(self):
        """Clustering molecules into groups based on molecular types."""
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