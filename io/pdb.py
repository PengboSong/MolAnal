# coding=utf-8

from gmx.io.baseio import GMXDomains, gmx_readline

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
