# coding=utf-8

import re

from gmx.io.baseio import GMXDomains, gmx_readline_sep


CSV_DOMAINS = {
    # Short Name : Long Name, loc, type, type name, aligned
    GMXDomains.ATOMID:  ("New atom serial number", 1, int, "interger", "right"),
    GMXDomains.OATOMID: ("Old atom serial number", 2, int, "interger", "right"),
    GMXDomains.ATOMNM:  ("Atom name", 3, str, "character", "left"),
    GMXDomains.RESNM:   ("Residual name", 4, str, "character", "right"),
    GMXDomains.RESID:   ("Residue sequence number", 5, int, "interger", "right"),
    GMXDomains.OCCUP:   ("Occupancy", 6, float, "real", "right"),
    GMXDomains.TEMP:    ("Temperature factor", 7, float, "real", "right"),
    GMXDomains.ELEMENT: ("Element symbol", 8, str, "character", "right"),
}


def readline_csv(line, row):
    return gmx_readline_sep(line, row, domains=CSV_DOMAINS)
