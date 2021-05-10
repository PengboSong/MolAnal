# coding=utf-8

import os

from gmx.io.baseio import GMXDomains
from gmx.io.csv import readline_csv
from gmx.io.pdb import readline_pdb
from gmx.structure.atom_matrix import AtomMatrix


def pdb2gro(pdb_path, csv_path, gro_path, maxsize=209715200):
    # Initialize dict
    data = {}
    modif = {}

    # Load pdb file
    if not os.path.isfile(pdb_path) or os.path.splitext(pdb_path)[-1] != ".pdb":
        raise IOError(
            "Can not load file {}. Expect a PDB file.".format(pdb_path))
    with open(pdb_path, 'r') as f:
        row = 0
        for line in f.readlines(maxsize):
            row += 1
            if line[0:6] in ("ATOM  ", "HETATM"):
                line_data = readline_pdb(line, row)
                old_atomid = line_data.get(GMXDomains.ATOMID)
                data.update({old_atomid: line_data})

    # Load csv file recording completion lines
    if not os.path.isfile(csv_path) or os.path.splitext(csv_path)[-1] != ".csv":
        raise IOError(
            "Can not load file {}. Expect a CSV file.".format(csv_path))
    with open(csv_path, 'r') as f:
        row = 0
        for line in f.readlines(maxsize):
            row += 1
            line_modif = readline_csv(line, row)
            new_atomid = line_modif.get(GMXDomains.ATOMID)
            modif.update({new_atomid: line_modif})

    outmat = AtomMatrix()
    for newid in modif:
        newitem = modif[newid]
        olditem = data[newitem[GMXDomains.ATOMID]]
        outmat.append(atomnm=olditem[GMXDomains.ATOMNM],
                      element=olditem[GMXDomains.ELEMENT],
                      molid=newitem[GMXDomains.RESID],
                      molnm=newitem[GMXDomains.RESNM],
                      *olditem[GMXDomains.XYZ])
    outmat.to_gro(gro_path)
