import os
from numpy import array as numpyArray
from numpy import vstack as verticalStack
from gmx.io.pdb import readline_pdb
from gmx.io.gro import readline_gro
from gmx.math.common import boxpara

def pdb2molmatrix(pdbpath, maxsize = 209715200):
	# Initialize molecular matrix dict
	mols = {}
	# Read input file
	with open(pdbpath, 'r') as f:
		# Start counting number of lines
		rown = 0
		former = (0,"")
		# pdb file
		if os.path.splitext(pdbpath)[1] == ".pdb":
			title = ""
			atomn = 0
			for line in f.readlines(maxsize):
				if line[0:6] in ("ATOM  ", "HETATM"):
					line_data = readline_pdb(line, rown)
					# First molecule
					if not former[0] and not former[1]:
						mol_id = line_data.get("residue_sequence_number")
						mol_name = line_data.get("residual_name")
						mol_matrix = {"name":mol_name}
						atom, coordinate, velocity = [], [], []
					# Start a new molecule
					elif former[0] != line_data.get("residue_sequence_number") or former[1] != line_data.get("residual_name"):
						mol_matrix.update({"atom":atom, "coordinate":numpyArray(coordinate), "velocity":numpyArray(velocity)})
						mols.update({mol_id:mol_matrix})
						mol_id = line_data.get("residue_sequence_number")
						mol_name = line_data.get("residual_name")
						mol_matrix = {"name":mol_name}
						atom, coordinate, velocity = [], [], []
					# Common lines
					atom.append(line_data.get("atom_name"))
					coordinate.append([line_data.get("coordinate_x")/10, line_data.get("coordinate_y")/10, line_data.get("coordinate_z")/10])
					# No velocity data in pdb file
					velocity = [[0., 0., 0.]]
					former = (line_data.get("residue_sequence_number"), line_data.get("residual_name"))
			# Complete the last molecular matrix
			mol_matrix.update({"atom":atom, "coordinate":numpyArray(coordinate), "velocity":numpyArray(velocity)})
			mols.update({mol_id:mol_matrix})
		# gro file
		elif os.path.splitext(pdbpath)[1] == ".gro":
			# Start counting and read first two lines
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
					# First molecule
					if not former[0] and not former[1]:
						mol_id = line_data.get("mol_id")
						mol_name = line_data.get("mol_name")
						mol_matrix = {"name":mol_name}
						atom, coordinate, velocity = [], [], []
					# Start a new molecule
					elif former[0] != line_data.get("mol_id") or former[1] != line_data.get("mol_name"):
						mol_matrix.update({"atom":atom, "coordinate":numpyArray(coordinate), "velocity":numpyArray(velocity)})
						mols.update({mol_id:mol_matrix})
						mol_id = line_data.get("mol_id")
						mol_name = line_data.get("mol_name")
						mol_matrix = {"name":mol_name}
						atom, coordinate, velocity = [], [], []
					# Common lines
					atom.append(line_data.get("atom_name"))
					coordinate.append([line_data.get("x"), line_data.get("y"), line_data.get("z")])
					velocity.append([line_data.get("vx"), line_data.get("vy"), line_data.get("vz")])
					former = (line_data.get("mol_id"), line_data.get("mol_name"))
			# Complete the last molecular matrix
			mol_matrix.update({"atom":atom, "coordinate":numpyArray(coordinate), "velocity":numpyArray(velocity)})
			mols.update({mol_id:mol_matrix})
		else:
			raise ValueError("Can not load file that does not end with \"pdb\" or \"gro\".")
	return mols

def molmatrix2pdb(mols, pdbpath, writevelo = False):
	# Start counting
	atomn, moln = 0, 0
	writeAtoms = []
	# Get the solvate box size
	xbox, ybox, zbox = boxpara(verticalStack([mol.get("coordinate") for mol in mols.values()]))
	# Write gro file
	if os.path.splitext(pdbpath)[1] == ".gro":
		for j in range(1, len(mols)+1):
			moln += 1
			mol = mols.get(j)
			molid = mol.get("id")
			moltype = mol.get("name")
			atom = mol.get("atom")
			coordinate = mol.get("coordinate")
			velocity = mol.get("velocity")
			for i in range(len(atom)):
				atomn += 1
				# Write velocity
				if writevelo is True:
					writeAtoms.append('{0:>5}'.format(moln) + '{0:<4}'.format(moltype) + " "*2 +'{0:>4}'.format(atom[i]) + '{0:>5}'.format(atomn) + '{0:>8.3f}'.format(coordinate[i][0]) + '{0:>8.3f}'.format(coordinate[i][1]) + '{0:>8.3f}'.format(coordinate[i][2]) + '{0:>8.3f}'.format(velocity[i][0]) + '{0:>8.3f}'.format(velocity[i][1]) + '{0:>8.3f}'.format(velocity[i][2]))
				# No velocity information
				else:
					writeAtoms.append('{0:>5}'.format(moln) + '{0:<4}'.format(moltype) + " "*2 +'{0:>4}'.format(atom[i]) + '{0:>5}'.format(atomn) + '{0:>8.3f}'.format(coordinate[i][0]) + '{0:>8.3f}'.format(coordinate[i][1]) + '{0:>8.3f}'.format(coordinate[i][2]))
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
		for j in range(1, len(mols)+1):
			moln += 1
			mol = mols.get(j)
			molid = mol.get("id")
			moltype = mol.get("name")
			atom = mol.get("atom")
			coordinate = mol.get("coordinate")
			for i in range(len(atom)):
				atomn += 1
				writeAtoms.append('{0:<6}'.format("HETATM") + '{0:>5}'.format(atomn) + " " + '{0:<4}'.format(atom[i]) + " " + '{0:>3}'.format(moltype) + " "*2 + '{0:>4}'.format(moln) + " "*4 + '{0:>8.3f}'.format(coordinate[i][0]*10) + '{0:>8.3f}'.format(coordinate[i][1]*10) + '{0:>8.3f}'.format(coordinate[i][2]*10) + '{0:>6.2f}'.format(1) + '{0:>6.2f}'.format(0) + " "*10 + '{0:>3}'.format(atom[i][0]))
			with open(pdbpath, "w") as f:
				f.writelines(line + '\n' for line in writeAtoms)
	else:
		raise ValueError("Can not write file that does not end with \"pdb\" or \"gro\".")
