# coding=utf-8

import os.path
import numpy as np

from gmx.io.baseio import GMXDomains, gmx_readline
from gmx.io.csv import readline_csv
from gmx.other.input_func import convert_input_type
from gmx.math.common import boxpara
from gmx.structure.atom_matrix import AtomMatrix

PDB_DOMAINS = {
    # Short Name : Long Name, Start loc, End loc, type, type name, aligned
    GMXDomains.RECORD:  ("Record type", 0, 6, str, "character", "left"),
    GMXDomains.ATOMID:  ("Atom serial number", 6, 11, int, "interger", "right"),
    GMXDomains.ATOMNM:  ("Atom name", 12, 16, str, "character", "left"),
    GMXDomains.LOC:     ("Alternate location indicator", 16, 17, str, "character", "none"),
    GMXDomains.RESNM:   ("Residual name", 17, 20, str, "character", "right"),
    GMXDomains.CHAIN:   ("Chain identifier", 21, 22, str, "character", "none"),
    GMXDomains.RESID:   ("Residue sequence number", 22, 26, int, "interger", "right"),
    GMXDomains.CODE:    ("Code for insertions of residuals", 26, 27, str, "character", "none"),
    GMXDomains.X:       ("Coordinate X(A)", 30, 38, float, "real", "right"),
    GMXDomains.Y:       ("Coordinate Y(A)", 38, 46, float, "real", "right"),
    GMXDomains.Z:       ("Coordinate Z(A)", 46, 54, float, "real", "right"),
    GMXDomains.OCCUP:   ("Occupancy", 54, 60, float, "real", "right"),
    GMXDomains.TEMP:    ("Temperature factor", 60, 66, float, "real", "right"),
    GMXDomains.SEGMENT: ("Segment identifier", 72, 76, str, "character", "left"),
    GMXDomains.ELEMENT: ("Element symbol", 76, 78, str, "character", "right"),
}


def readline_pdb(line, row):
    line = gmx_readline(line, row, domains=PDB_DOMAINS)
    # Convert x, y, z from angstorm to nanometer
    line[GMXDomains.X] *= .1
    line[GMXDomains.Y] *= .1
    line[GMXDomains.Z] *= .1
    # Pack coordinate XYZ
    line[GMXDomains.XYZ] = [line[GMXDomains.X],
                            line[GMXDomains.Y],
                            line[GMXDomains.Z]]
    return line


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
