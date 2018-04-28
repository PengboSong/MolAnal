import os
from numpy import array as numpyArray
from numpy import vstack as verticalStack
from gmx.io.pdb import readline_pdb
from gmx.io.gro import readline_gro
from gmx.math.common import boxpara

def pdb2atomatrix(pdbpath, maxsize = 209715200):
	# Initialize atom matrix dict
	atoms = {}
	# Read input file
	with open(pdbpath, "r") as f:
		# Start counting number of lines
		rown = 0
		# pdb file
		if os.path.splitext(pdbpath)[1] == ".pdb":
			for line in f.readlines(maxsize):
				if line[0:6] in ("ATOM  ", "HETATM"):
					line_data = readline_pdb(line, rown)
					atom_matrix = {"name":line_data.get("atom_name"), "resid":line_data.get("residue_sequence_number"), "resname":line_data.get("residual_name"), "coordinate":numpyArray([line_data.get("coordinate_x")/10, line_data.get("coordinate_y")/10, line_data.get("coordinate_z")/10]), "velocity":numpyArray([0., 0., 0.])}
					atoms.update({line_data.get("atom_serial_number"):atom_matrix})
		# gro file
		elif os.path.splitext(pdbpath)[1] == ".gro":
			rown += 1
			title = f.readline().strip()
			rown += 1
			atomn = f.readline().strip()
			# Second line should be an integer equals to total atom numbers
			try:
				atomn = int(atomn)
			except Exception as e:
				raise ValueError("Incorrect gromacs file format. Please check file at path %s." % pdbpath)
			for line in f.readlines(maxsize):
				rown += 1
				# Line 'atomn + 3' should be solvate box parameters
				if rown < atomn + 3:
					line_data = readline_gro(line, rown)
					atom_matrix = {"name":line_data.get("atom_name"), "resid":line_data.get("mol_id"), "resname":line_data.get("mol_name"), "coordinate":numpyArray([line_data.get("x"), line_data.get("y"), line_data.get("z")]), "velocity":numpyArray([line_data.get("vx"), line_data.get("vy"), line_data.get("vz")])}
					atoms.update({line_data.get("atom_id"):atom_matrix})
		else:
			raise ValueError("Can not load file that does not end with \"pdb\" or \"gro\".")
	return atoms


def atomatrix2pdb(atoms, pdbpath, writevelo = False):
	# Start counting
	atomn, moln = 0, 0
	writeAtoms = []
	# Get the solvate box size
	xbox, ybox, zbox = boxpara(verticalStack([atom.get("coordinate") for atom in self.atoms.values()]))
	# Write gro file
	if os.path.splitext(pdbpath)[1] == ".gro":
		for j in range(1, len(atoms)+1):
			atomn += 1
			atom = atoms.get(j)
			atomtype = atom.get("name")
			molid = atom.get("resid")
			moltype = atom.get("resname")
			coordinate = atom.get("coordinate")
			velocity = atom.get("velocity")
			# New molecule
			if j == 1 or (j > 1 and (moltype != atoms.get(j-1).get("resname") or molid != atoms.get(j-1).get("resid"))):
				moln += 1
			if writevelo is True:
				# Need completion
				writeAtoms.append()
			else:
				# Need completion
				writeAtoms.append()
		with open(pdbpath, "w") as f:
			# Default title
			f.write("GROtesk MACabre and Sinister" + '\n')
			# Atom numbers
			f.write('{0:>5}'.format(atomn) + '\n')
			f.writelines(line + '\n' for line in writeAtoms)
			# Solvate box parameters
			f.write('{0:>10.5f}'.format(xbox) + '{0:>10.5f}'.format(ybox) + '{0:>10.5f}'.format(zbox) + '\n')
	# Write pdb file
	elif os.path.splitext(pdbpath)[1] == ".pdb":
		# Need completion
		pass
	else:
		raise ValueError("Can not write file that does not end with \"pdb\" or \"gro\".")
