# coding=utf-8

from gmx.io.baseio import GMXDomains, gmx_readline


GRO_DOMAINS = {
    # Short Name : Long Name, Start loc, End loc, type, type name, aligned
    GMXDomains.MOLID:   ("Molecular id", 0, 5, int, "interger", "right"),
    GMXDomains.MOLNM:   ("Molecular name", 5, 8, str, "character", "right"),
    GMXDomains.CHAIN:   ("Chain identifier", 9, 10, str, "character", "none"),
    GMXDomains.ATOMNM:  ("Atom name", 11, 15, str, "character", "right"),
    GMXDomains.ATOMID:  ("Atom id", 15, 20, int, "interger", "right"),
    GMXDomains.X:       ("Coordinate X(nm)", 20, 28, float, "real", "right"),
    GMXDomains.Y:       ("Coordinate Y(nm)", 28, 36, float, "real", "right"),
    GMXDomains.Z:       ("Coordinate Z(nm)", 36, 44, float, "real", "right"),
    GMXDomains.VX:      ("Velocity X", 44, 52, float, "real", "right"),
    GMXDomains.VY:      ("Velocity Y", 52, 60, float, "real", "right"),
    GMXDomains.VZ:      ("Velocity Z", 60, 68, float, "real", "right"),
}


def readline_gro(line, row):
    line = gmx_readline(line, row, domains=GRO_DOMAINS)
    # Pack coordinate and velocity
    line[GMXDomains.XYZ] = [line[GMXDomains.X],
                            line[GMXDomains.Y],
                            line[GMXDomains.Z]]
    if GMXDomains.VX in line and GMXDomains.VY in line and GMXDomains.VZ in line:
        line[GMXDomains.VXYZ] = [line[GMXDomains.VX],
                                line[GMXDomains.VY],
                                line[GMXDomains.VZ]]
    return line
