import os.path
import numpy as np

from gmx.io.csv import readline_csv_spec
from gmx.other.input_func import convert_input_type
from gmx.math.common import boxpara

def readline_pdb(line, row):
	assert isinstance(row, int), "Invalid row number."
	assert isinstance(line, str), "Can not parse the given line at row %d." % row
	# Read by column number

	# Record type: 1-6, character, left
	record = convert_input_type(line[0:6].strip())
	assert isinstance(record, str), "Wrong data type. Record type at row %d column 1-6 is not characters." % row

	# Atom serial number: 7-11, integer, right
	atom_id = convert_input_type(line[6:11].strip())
	assert isinstance(atom_id, int), "Wrong data type. Atom serial number at row %d column 7-11 is not an integer." % row

	# Atom name: 13-16, character, left
	atom_name = convert_input_type(line[12:16].strip())
	assert isinstance(atom_name, str), "Wrong data type. Atom name at row %d column 13-16 is not characters." % row

	# Alternate location indicator: 17, character
	location = convert_input_type(line[16:17].strip())
	if location:
		assert isinstance(location, str), "Wrong data type. Alternate location at row %d column 17 is not a character." % row

	# Residual name: 18-20, character, right
	res_name = convert_input_type(line[17:20].strip())
	if res_name:
		assert isinstance(res_name, str), "Wrong data type. Residual name at row %d column 18-20 is not characters." % row

	# Chain identifier: 22, character
	chain = convert_input_type(line[21:22].strip())
	if chain:
		assert isinstance(chain, str), "Wrong data type. Chain identifier at row %d column 22 is not characters." % row

	# Residue sequence number: 23-26, interger, right
	res_id = convert_input_type(line[22:26].strip())
	if res_id:
		assert isinstance(res_id, int), "Wrong data type. Residue sequence number at row %d column 23-26 is not an integer." % row

	# Code for insertions of residuals: 27, character
	code = convert_input_type(line[26:27].strip())
	if code:
		assert isinstance(code, str), "Wrong data type. Code for insertions of residuals at row %d column 27 is not a character." % row

	# Coordinate X(A): 31-38, real(8:3), right
	x = convert_input_type(line[30:38].strip())
	assert isinstance(x, float), "Wrong data type. Coordinate X at row %d column 31-38 is not a real number." % row
	# Coordinate Y(A): 39-46, real(8:3), right
	y = convert_input_type(line[38:46].strip())
	assert isinstance(y, float), "Wrong data type. Coordinate Y at row %d column 39-46 is not a real number." % row
	# Coordinate Z(A): 47-54, real(8:3), right
	z = convert_input_type(line[46:54].strip())
	assert isinstance(z, float), "Wrong data type. Coordinate Z at row %d column 47-54 is not a real number." % row

	# Occupancy: 55-60, real(6:2), right
	occup = convert_input_type(line[54:60].strip())
	if occup:
		assert isinstance(occup, float), "Wrong data type. Occupancy at row %d column 55-60 is not a real number." % row

	# Temperature factor: 61-66, real(6:2), right
	temp_factor = convert_input_type(line[60:66].strip())
	if temp_factor:
		assert isinstance(temp_factor, float), "Wrong data type. Temperature factor at row %d column 61-66 is not a real number." % row

	# Segment identifier: 73-76, character, left
	segment = convert_input_type(line[72:76].strip())
	if segment:
		assert isinstance(segment, str), "Wrong data type. Segment identifier at row %d column 73-76 is not characters." % row

	# Element symbol: 77-78, character, right
	element = convert_input_type(line[76:78].strip())
	if element:
		assert isinstance(element, str), "Wrong data type. Element symbol at row %d column 77-78 is not characters." % row

	return {
		"record_type":record,
		"atom_serial_number":atom_id,
		"atom_name":atom_name,
		"alternate_location_identifier":location,
		"residual_name":res_name,
		"chain_identifier":chain,
		"residue_sequence_number":res_id,
		"code_for_insertions_of_residues":code,
		# Convert coordinate from angstrom to nanometer
		"coordinate":np.array([x, y, z])/10,
		"occupancy":occup,
		"temperature_factor":temp_factor,
		"segment_identifier":segment,
		"element_symbol":element
	}

def pdb2atom_matrix(pdb_path, maxsize = 209715200):
	# Initialize atom matrix dict
	atoms = {}
	# Load pdb file
	if os.path.isfile(pdb_path) and os.path.splitext(pdb_path)[1] == ".pdb":
		with open(pdb_path, "r") as f:
			# Start counting number of lines
			rown = 0
			for line in f.readlines(maxsize):
				if line[0:6] in ("ATOM  ", "HETATM"):
					line_data = readline_pdb(line, rown)
					atom_matrix = {
						"name":line_data.get("atom_name"),
						"resid":line_data.get("residue_sequence_number"),
						"resname":line_data.get("residual_name"),
						"coordinate":line_data.get("coordinate"),
						"velocity":np.array([0., 0., 0.])
					}
					atoms.update({line_data.get("atom_serial_number"):atom_matrix})
	else:
		raise IOError("Can not load file %s." % pdb_path)
	return atoms

def pdb2coord_matrix(pdb_path, maxsize = 209715200):
	# Initialize coordinate matrix
	coords = []
	# Load pdb file
	if os.path.isfile(pdb_path) and os.path.splitext(pdb_path)[1] == ".pdb":
		with open(pdb_path, 'r') as f:
			rown = 0
			for line in f.readlines(maxsize):
				rown += 1
				if line[0:6] in ("ATOM  ", "HETATM"):
					coords.append(readline_pdb(line, rown).get("coordinate"))
	else:
		raise IOError("Can not load file %s." % pdb_path)
	coords = np.vstack(coords)
	return coords

def pdb2mol_matrix(pdb_path, maxsize = 209715200):
	# Initialize molecular matrix dict
	mols = {}
	# Load pdb file
	if os.path.isfile(pdb_path) and os.path.splitext(pdb_path)[1] == ".pdb":
		with open(pdb_path, 'r') as f:
			# Start counting number of lines
			rown, atomn = 0, 0
			former = (0,"")
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
						mol_matrix.update({"atom":atom, "coordinate":np.array(coordinate), "velocity":np.array(velocity)})
						mols.update({mol_id:mol_matrix})
						mol_id = line_data.get("residue_sequence_number")
						mol_name = line_data.get("residual_name")
						mol_matrix = {"name":mol_name}
						atom, coordinate, velocity = [], [], []
					# Common lines
					atom.append(line_data.get("atom_name"))
					coordinate.append(line_data.get("coordinate"))
					# No velocity data in pdb file
					velocity.append([0., 0., 0.])
					former = (line_data.get("residue_sequence_number"), line_data.get("residual_name"))
			# Complete the last molecular matrix
			mol_matrix.update({"atom":atom, "coordinate":np.array(coordinate), "velocity":np.array(velocity)})
			mols.update({mol_id:mol_matrix})
	else:
		raise IOError("Can not load file %s." % pdb_path)
	return mols

def pdb2gro(pdb_path, csv_path, gro_path, maxsize = 209715200):
	# Initialize dict
	data = {}
	modif = {}

	# Load pdb file
	if os.path.isfile(pdb_path) and os.path.splitext(pdb_path)[1] == ".pdb":
		with open(pdb_path, 'r') as f:
			row = 0
			for line in f.readlines(maxsize):
				row += 1
				if line[0:6] in ("ATOM  ", "HETATM"):
					line_data = readline_pdb(line, row)
					data.update({line_data.get("atom_serial_number"):line_data})
	else:
		raise IOError("Can not load file %s." % pdb_path)

	# Load csv file recording completion lines
	if os.path.isfile(csv_path) and os.path.splitext(csv_path)[1] == ".csv":
		with open(csv_path, 'r') as f:
			row = 0
			for line in f.readlines(maxsize):
				row += 1
				line_modif = readline_csv_spec(line, row)
				if line_modif:
					modif.update({line_modif.get("new_atom_serial_number"):line_modif})
	else:
		raise IOError("Can not load file %s. Expect a .csv file." % csv_path)

	# Write parts
	lxyz, lout = [], []
	for newid in modif.keys():
		modif_term = modif.get(newid)
		out_term = data.get(modif_term.get("atom_serial_number")).copy()
		out_term.update(modif_term)
		out_term.update({"molecular_id":modif_term.get("residue_sequence_number")})
		lout.append(out_term)
		lxyz.append(out_term.get("coordinate"))

	# Get the solvate box parameters
	xbox, ybox, zbox = boxpara(np.array(lxyz))
	with open(gro_path, 'w') as f:
		# Default title
		f.write("GROtesk MACabre and Sinister" + '\n')
		# Total atom numbers
		f.write('{0:>5}'.format(len(lout)) + '\n')
		new_lines = []
		for term in lout:
			this_x, this_y, this_z = term.get("coordinate")
			term_line = '{0:>5}'.format(term.get("molecular_id") or " ")
			term_line += '{0:<4}'.format(term.get("residual_name") or " ")
			term_line += '{0:1}'.format(term.get("chain_identifier") or " ")
			term_line += '{0:>4}'.format(term.get("atom_name") or " ")
			term_line += '{0:>5}'.format(term.get("atom_serial_number") or " ")
			term_line += '{0:>8.3f}'.format(this_x or " ") + '{0:>8.3f}'.format(this_y or " ") + '{0:>8.3f}'.format(this_z or " ")
			new_lines.append(term_line)
		f.writelines(line + '\n' for line in term_lines)
		# Solvate box parameters
		f.write('{0:>10.5f}'.format(xbox) + '{0:>10.5f}'.format(ybox) + '{0:>10.5f}'.format(zbox) + '\n')

def mol_matrix2pdb(mols, file_name = "frame"):
	# Start counting
	atomn, moln = 0, 0
	write_atoms = []
	# Write pdb file
	for j in range(len(mols)):
		moln += 1
		mol = mols.get(j + 1)
		molid = mol.get("id")
		moltype = mol.get("name")
		atom = mol.get("atom")
		coordinate = mol.get("coordinate")
		for i in range(len(atom)):
			x, y, z = coordinate[i] * 10
			atomn += 1
			atom_line = '{0:<6}'.format("HETATM")
			atom_line += '{0:>5}'.format(atomn)
			atom_line += ' '
			atom_line += '{0:<4}'.format(atom[i])
			atom_line += ' '
			atom_line += '{0:>3}'.format(moltype)
			atom_line += ' ' * 2
			atom_line += '{0:>4}'.format(moln)
			atom_line += ' ' * 4
			atom_line += '{0:>8.3f}'.format(x)
			atom_line += '{0:>8.3f}'.format(y)
			atom_line += '{0:>8.3f}'.format(z)
			atom_line += '{0:>6.2f}'.format(1)
			atom_line += '{0:>6.2f}'.format(20)
			atom_line += ' ' * 10
			atom_line += '{0:>3}'.format(atom[i][0])
			write_atoms.append(atom_line)
	with open(file_name + ".pdb", "w") as f:
		f.writelines(line + '\n' for line in write_atoms)
