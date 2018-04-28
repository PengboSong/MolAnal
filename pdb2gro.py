import os
from numpy import array as numpyArray
from gmx.io.csv import readline_csv_spec
from gmx.io.pdb import readline_pdb
from gmx.math.common import boxpara

def pdb2gro(pdbpath, csvpath, writepath, maxsize = 209715200):
    # Initialize dict
    data = {}
    modif = {}
    # Load pdb file
    with open(pdbpath, 'r') as f:
        row = 0
        for line in f.readlines(maxsize):
            row += 1
            if line[0:6] in ("ATOM  ", "HETATM"):
                line_data = readline_pdb(line, row)
                data.update({line_data.get("atom_serial_number"):line_data})
    # Load csv file recording completion lines
    with open(csvpath, 'r') as f:
        row = 0
        for line in f.readlines(maxsize):
            row += 1
            line_modif = readline_csv_spec(line, row)
            if line_modif:
                modif.update({line_modif.get("new_atom_serial_number"):line_modif})
    # Write parts
    lxyz, lout = [], []
    for newid in modif.keys():
        modifterm = modif.get(newid)
        outerm = data.get(modifterm.get("atom_serial_number")).copy()
        outerm.update(modifterm)
        outerm.update({"molecular_id":modifterm.get("residue_sequence_number")})
        lout.append(outerm)
        lxyz.append((outerm.get("coordinate_x"), outerm.get("coordinate_y"), outerm.get("coordinate_z")))
    # Get the
    xbox, ybox, zbox = boxpara(numpyArray(lxyz))
    with open(writepath, 'w') as f:
        if os.path.splitext(pdbpath)[1] == ".pdb":
            # Need Completion
            pass
        if os.path.splitext(pdbpath)[1] == ".gro":
            # Default title
            f.write("GROtesk MACabre and Sinister" + '\n')
            # Total atom numbers
            f.write('{0:>5}'.format(len(lout)) + '\n')
            for term in lout:
                f.write('{0:>5}'.format(term.get("molecular_id") or " ") + '{0:<4}'.format(term.get("residual_name") or " ") + '{0:1}'.format(term.get("chain_identifier") or " ") + '{0:>4}'.format(term.get("atom_name") or " ") + '{0:>5}'.format(term.get("atom_serial_number") or " ") + '{0:>8}'.format(format(term.get("coordinate_x")/10, '.3f') or " ") + '{0:>8}'.format(format(term.get("coordinate_y")/10, '.3f') or " ") + '{0:>8}'.format(format(term.get("coordinate_z")/10, '.3f') or " ") + '\n')
            # Solvate box parameters
            f.write('{0:>10.5f}'.format(xbox/10) + '{0:>10.5f}'.format(ybox/10) + '{0:>10.5f}'.format(zbox/10) + '\n')
