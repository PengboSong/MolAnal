# Import frequent used modules
# os: file reading and writeing
import os
# math: common mathematics operations
import math
# numpy: matrix and scientific calculations
import numpy as np
# Functions extracted from older edition
def check_input_type(inp, typ, restrain = None):
	if not isinstance(inp, typ):
		return False
	if restrain:
		if isinstance(restrain, str):
			exec("result = " + repr(inp) + restrain)
			return result
		if isinstance(restrain, (tuple, list)):
			results = []
			for term in restrain:
				exec("results.append(" + repr(inp) + term + ")")
			return all(results)
	else:
		return True
def verify_type(typ, desc, restrain = None, msg = "Invalid input."):
	val = convert_input_type(input(desc))
	while not check_input_type(val, typ, restrain):
		print(msg)
		val = convert_input_type(input(desc))
	return val
def verify_selection(ls, desc, msg = "Invalid option."):
	if not isinstance(ls, (list, dict, tuple)):
		raise ValueError("Bad data type. Expect a list, a dict or a tuple.")
	if isinstance(ls, dict):
		for i in range(len(ls)):
			print(format(i+1, 'd') + ':' + (ls.get(i+1) or "None"))
	inp = convert_input_type(input(desc))
	while inp not in ls:
		print(msg)
		inp = convert_input_type(input(desc))
	return inp
def convert_input_type(string):
	# Only handle strings
	if isinstance(string, str):
		# First, judge whether the string is an integer
		try:
			integer = int(string)
			return integer
		except ValueError as e:
			# Second, judge whether it is a floating number
			try:
				floating = float(string)
				return floating
			except ValueError as e:
				# If it is not a number, regard it as a string
				try:
					# Caution: No check
					eval(string)
				except Exception:
					return string
				else:
					return eval(string)
	else:
		return string
def readline_csv(line, split_symbol = ';', comment_symbol = '#'):
	if not line.startswith(comment_symbol):
		return line.split(comment_symbol, 1)[0].strip().split(split_symbol)
	else:
		return None
def readline_csv_spec(line, row):
	# Newline; Oldline; AtomName; ResType; ResId; Occup; TempFactor; SegName
	if readline_csv(line):
		new_id, old_id, atom_name, res_type, res_id, occup, temp_factor, segment = readline_csv(line)
		if isinstance(convert_input_type(new_id), int):
			new_id = convert_input_type(new_id)
		else:
			raise ValueError("Wrong data type. Term Newline at row " + format(row, 'd') + " column 1 is not an integer.")
		if new_id != 0:
			if isinstance(convert_input_type(old_id), int):
				old_id = convert_input_type(old_id)
			else:
				raise ValueError("Wrong data type. Term Oldline at row " + format(row, 'd') + " column 2 is not an integer.")
			if isinstance(convert_input_type(atom_name), str):
				atom_name = convert_input_type(atom_name)
			else:
				raise ValueError("Wrong data type. Term AtomName at row " + format(row, 'd') + " column 3 is not characters.")
			if isinstance(convert_input_type(res_type), str):
				res_type = convert_input_type(res_type)
			else:
				raise ValueError("Wrong data type. Term ResType at row " + format(row, 'd') + " column 4 is not characters.")
			if isinstance(convert_input_type(res_id), int):
				res_id = convert_input_type(res_id)
			else:
				raise ValueError("Wrong data type. Term ResId at row " + format(row, 'd') + " column 5 is not an integer.")
			if isinstance(convert_input_type(occup), float):
				occup = convert_input_type(occup)
			else:
				raise ValueError("Wrong data type. Term Occup at row " + format(row, 'd') + " column 6 is not a real number.")
			if isinstance(convert_input_type(temp_factor), float):
				temp_factor = convert_input_type(temp_factor)
			else:
				raise ValueError("Wrong data type. Term TempFactor at row " + format(row, 'd') + " column 7 is not a real number.")
			if isinstance(convert_input_type(segment), str):
				segment = convert_input_type(segment)
			else:
				raise ValueError("Wrong data type. Term SegName at row " + format(row, 'd') + " column 8 is not characters.")
			return {"new_atom_serial_number":new_id, "atom_serial_number":old_id, "atom_name":atom_name, "residual_name":res_type, "residue_sequence_number":res_id, "occupancy":occup, "temperature_factor":temp_factor, "element_symbol":segment}
		else:
			return None
	else:
		return None
def readline_pdb(line, row):
	# Read by column number
	# Record type: 1-6, character, left
	record = convert_input_type(line[0:6].strip())
	if not isinstance(record, str):
		raise ValueError("Wrong data type. Record type at row " + format(row, 'd') + " column 1-6 is not characters.")
	# Atom serial number: 7-11, integer, right
	atom_id = convert_input_type(line[6:11].strip())
	if not isinstance(atom_id, int):
		raise ValueError("Wrong data type. Atom serial number at row " + format(row, 'd') + " column 7-11 is not an integer.")
	# Atom name: 13-16, character, left
	atom_name = convert_input_type(line[12:16].strip())
	if not isinstance(atom_name, str):
		raise ValueError("Wrong data type. Atom name at row " + format(row, 'd') + " column 13-16 is not characters.")
	# Alternate location indicator: 17, character
	location = convert_input_type(line[16:17].strip())
	if location and not isinstance(location, str):
		raise ValueError("Wrong data type. Alternate location at row " + format(row, 'd') + " column 17 is not a character.")
	# Residual name: 18-20, character, right
	res_name = convert_input_type(line[17:20].strip())
	if res_name and not isinstance(res_name, str):
		raise ValueError("Wrong data type. Residual name at row " + format(row, 'd') + " column 18-20 is not characters.")
	# Chain identifier: 22, character
	chain = convert_input_type(line[21:22].strip())
	if chain and not isinstance(chain, str):
		raise ValueError("Wrong data type. Chain identifier at row " + format(row, 'd') + " column 22 is not characters.")
	# Residue sequence number: 23-26, interger, right
	res_id = convert_input_type(line[22:26].strip())
	if res_id and not isinstance(res_id, int):
		raise ValueError("Wrong data type. Residue sequence number at row " + format(row, 'd') + " column 23-26 is not an integer.")
	# Code for insertions of residuals: 27, character
	code = convert_input_type(line[26:27].strip())
	if code and not isinstance(code, str):
		raise ValueError("Wrong data type. Code for insertions of residuals at row " + format(row, 'd') + " column 27 is not a character.")
	# Coordinate X(A): 31-38, real(8:3), right
	x = convert_input_type(line[30:38].strip())
	if not isinstance(x, float):
		raise ValueError("Wrong data type. Coordinate X at row " + format(row, 'd') + " column 31-38 is not a real number.")
	# Coordinate Y(A): 39-46, real(8:3), right
	y = convert_input_type(line[38:46].strip())
	if not isinstance(y, float):
		raise ValueError("Wrong data type. Coordinate Y at row " + format(row, 'd') + " column 39-46 is not a real number.")
	# Coordinate Z(A): 47-54, real(8:3), right
	z = convert_input_type(line[46:54].strip())
	if not isinstance(z, float):
		raise ValueError("Wrong data type. Coordinate Z at row " + format(row, 'd') + " column 47-54 is not a real number.")
	# Occupancy: 55-60, real(6:2), right
	occup = convert_input_type(line[54:60].strip())
	if occup and not isinstance(occup, float):
		raise ValueError("Wrong data type. Occupancy at row " + format(row, 'd') + " column 55-60 is not a real number.")
	# Temperature factor: 61-66, real(6:2), right
	temp_factor = convert_input_type(line[60:66].strip())
	if temp_factor and not isinstance(temp_factor, float):
		raise ValueError("Wrong data type. Temperature factor at row " + format(row, 'd') + " column 61-66 is not a real number.")
	# Segment identifier: 73-76, character, left
	segment = convert_input_type(line[72:76].strip())
	if segment and not isinstance(segment, str):
		raise ValueError("Wrong data type. Segment identifier at row " + format(row, 'd') + " column 73-76 is not characters.")
	# Element symbol: 77-78, character, right
	element = convert_input_type(line[76:78].strip())
	if element and not isinstance(element, str):
		raise ValueError("Wrong data type. Element symbol at row " + format(row, 'd') + " column 77-78 is not characters.")
	return {"record_type":record, "atom_serial_number":atom_id, "atom_name":atom_name, "alternate_location_identifier":location, "residual_name":res_name, "chain_identifier":chain, "residue_sequence_number":res_id, "code_for_insertions_of_residues":code, "coordinate_x":x, "coordinate_y":y, "coordinate_z":z, "occupancy":occup, "temperature_factor":temp_factor, "segment_identifier":segment, "element_symbol":element}
def readline_gro(line, row):
	# Read by column number
	# Molecular id: 1-5, character, right
	mol_id = convert_input_type(line[0:5].strip())
	if not isinstance(mol_id, int):
		raise ValueError("Wrong data type. Molecular id at row " + format(row, 'd') + " column 1-5 is not an integer.")
	# Molecular name: 6-8, character, left
	mol_name = convert_input_type(line[5:8].strip())
	if not isinstance(mol_name, str):
		raise ValueError("Wrong data type. Molecular name at row " + format(row, 'd') + " column 6-8 is not characters.")
	# Chain identifier: 10, character
	chain = convert_input_type(line[9:10].strip())
	if chain and not isinstance(chain, str):
		raise ValueError("Wrong data type. Chain identifier at row " + format(row, 'd') + " column 10 is not a character.")
	# Atom name: 12-15, character, right
	atom_name = convert_input_type(line[11:15].strip())
	if not isinstance(atom_name, str):
		raise ValueError("Wrong data type. Atom name at row " + format(row, 'd') + " column 12-15 is not characters.")
	# Atom id: 16-20, integer, right
	atom_id = convert_input_type(line[15:20].strip())
	if not isinstance(atom_id, int):
		raise ValueError("Wrong data type. Atom id at row " + format(row, 'd') + " column 16-20 is not an integer.")
	# Coordinate X(A): 21-28, real(8:3), right
	x = convert_input_type(line[20:28].strip())
	if not isinstance(x, float):
		raise ValueError("Wrong data type. Coordinate X at row " + format(row, 'd') + " column 21-28 is not a real number.")
	# Coordinate Y(A): 29-36, real(8:3), right
	y = convert_input_type(line[28:36].strip())
	if not isinstance(y, float):
		raise ValueError("Wrong data type. Coordinate Y at row " + format(row, 'd') + " column 29-36 is not a real number.")
	# Coordinate Z(A): 37-44, real(8:3), right
	z = convert_input_type(line[36:44].strip())
	if not isinstance(z, float):
		raise ValueError("Wrong data type. Coordinate Z at row " + format(row, 'd') + " column 37-44 is not a real number.")
	if len(line.strip()) > 44:
		# Velocity Vx: 45-52, read(8:4), right
		vx = convert_input_type(line[44:52].strip())
		if not isinstance(vx, float):
			raise ValueError("Wrong data type. Velocity Vx at row " + format(row, 'd') + " column 45-52 is not a real number.")
		# Velocity Vy: 53-60, read(8:4), right
		vy = convert_input_type(line[52:60].strip())
		if not isinstance(vy, float):
			raise ValueError("Wrong data type. Velocity Vy at row " + format(row, 'd') + " column 53-60 is not a real number.")
		# Velocity Vz: 61-68, read(8:4), right
		vz = convert_input_type(line[60:68].strip())
		if not isinstance(vz, float):
			raise ValueError("Wrong data type. Velocity Vz at row " + format(row, 'd') + " column 61-68 is not a real number.")
	else:
		vx, vy, vz = 0., 0., 0.
	return {"mol_id":mol_id, "mol_name":mol_name, "chain":chain, "atom_name":atom_name, "atom_id":atom_id, "x":x, "y":y, "z":z, "vx":vx, "vy":vy, "vz":vz}
# Attention: extend should be greater than 1
def boxpara(arrayxyz, extend = 1.10):
	if isinstance(arrayxyz, np.ndarray) and arrayxyz.ndim == 2 and arrayxyz.shape[1] == 3 and isinstance(extend, (int, float)) and extend > 1:
		xmin, ymin, zmin = np.min(arrayxyz, axis = 0)
		xmax, ymax, zmax = np.max(arrayxyz, axis = 0)
		return (xmax-xmin) * (1+extend), (ymax-ymin) * (1+extend), (zmax-zmin) * (1+extend)
	else:
		raise ValueError("Get wrong paramters. \"arrayxyz\" should be a two-dim numpy.ndarray with three columns and \"extend\" should be greater than 1.")
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
		xbox, ybox, zbox = boxpara(np.array(lxyz))
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
					atom_matrix = {"name":line_data.get("atom_name"), "resid":line_data.get("residue_sequence_number"), "resname":line_data.get("residual_name"), "coordinate":np.array([line_data.get("coordinate_x")/10, line_data.get("coordinate_y")/10, line_data.get("coordinate_z")/10]), "velocity":np.array([0., 0., 0.])}
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
					atom_matrix = {"name":line_data.get("atom_name"), "resid":line_data.get("mol_id"), "resname":line_data.get("mol_name"), "coordinate":np.array([line_data.get("x"), line_data.get("y"), line_data.get("z")]), "velocity":np.array([line_data.get("vx"), line_data.get("vy"), line_data.get("vz")])}
					atoms.update({line_data.get("atom_id"):atom_matrix})
		else:
			raise ValueError("Can not load file that does not end with \"pdb\" or \"gro\".")
	return atoms
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
						mol_matrix.update({"atom":atom, "coordinate":np.array(coordinate), "velocity":np.array(velocity)})
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
			mol_matrix.update({"atom":atom, "coordinate":np.array(coordinate), "velocity":np.array(velocity)})
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
						mol_matrix.update({"atom":atom, "coordinate":np.array(coordinate), "velocity":np.array(velocity)})
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
			mol_matrix.update({"atom":atom, "coordinate":np.array(coordinate), "velocity":np.array(velocity)})
			mols.update({mol_id:mol_matrix})
		else:
			raise ValueError("Can not load file that does not end with \"pdb\" or \"gro\".")
	return mols
def atomatrix2pdb(atoms, pdbpath, writevelo = False):
	# Start counting
	atomn, moln = 0, 0
	watoms = []
	# Get the solvate box size
	xbox, ybox, zbox = boxpara(np.vstack([atom.get("coordinate") for atom in self.atoms.values()]))
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
				watoms.append()
			else:
				# Need completion
				watoms.append()
		with open(pdbpath, "w") as f:
			# Default title
			f.write("GROtesk MACabre and Sinister" + '\n')
			# Atom numbers
			f.write('{0:>5}'.format(atomn) + '\n')
			f.writelines(line + '\n' for line in watoms)
			# Solvate box parameters
			f.write('{0:>10.5f}'.format(xbox) + '{0:>10.5f}'.format(ybox) + '{0:>10.5f}'.format(zbox) + '\n')
	# Write pdb file
	elif os.path.splitext(pdbpath)[1] == ".pdb":
		# Need completion
		pass
	else:
		raise ValueError("Can not write file that does not end with \"pdb\" or \"gro\".")
def molmatrix2pdb(mols, pdbpath, writevelo = False):
	# Start counting
	atomn, moln = 0, 0
	watoms = []
	# Get the solvate box size
	xbox, ybox, zbox = boxpara(np.vstack([mol.get("coordinate") for mol in mols.values()]))
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
					watoms.append('{0:>5}'.format(moln) + '{0:<4}'.format(moltype) + " "*2 +'{0:>4}'.format(atom[i]) + '{0:>5}'.format(atomn) + '{0:>8.3f}'.format(coordinate[i][0]) + '{0:>8.3f}'.format(coordinate[i][1]) + '{0:>8.3f}'.format(coordinate[i][2]) + '{0:>8.3f}'.format(velocity[i][0]) + '{0:>8.3f}'.format(velocity[i][1]) + '{0:>8.3f}'.format(velocity[i][2]))
				# No velocity information
				else:
					watoms.append('{0:>5}'.format(moln) + '{0:<4}'.format(moltype) + " "*2 +'{0:>4}'.format(atom[i]) + '{0:>5}'.format(atomn) + '{0:>8.3f}'.format(coordinate[i][0]) + '{0:>8.3f}'.format(coordinate[i][1]) + '{0:>8.3f}'.format(coordinate[i][2]))
		with open(pdbpath, "w") as f:
			# Default title
			f.write("GROtesk MACabre and Sinister" + '\n')
			# Atom numbers
			f.write('{0:>5}'.format(atomn) + '\n')
			f.writelines(line + '\n' for line in watoms)
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
				watoms.append('{0:<6}'.format("HETATM") + '{0:>5}'.format(atomn) + " " + '{0:<4}'.format(atom[i]) + " " + '{0:>3}'.format(moltype) + " "*2 + '{0:>4}'.format(moln) + " "*4 + '{0:>8.3f}'.format(coordinate[i][0]*10) + '{0:>8.3f}'.format(coordinate[i][1]*10) + '{0:>8.3f}'.format(coordinate[i][2]*10) + '{0:>6.2f}'.format(1) + '{0:>6.2f}'.format(0) + " "*10 + '{0:>3}'.format(atom[i][0]))
			with open(pdbpath, "w") as f:
				f.writelines(line + '\n' for line in watoms)
	else:
		raise ValueError("Can not write file that does not end with \"pdb\" or \"gro\".")
# Group molecules
def cluster(mols):
	groups = {}
	for molid in range(1, len(mols)+1):
		moltype = mols.get(molid).get("name")
		# New molecular group
		if not groups:
			groups.update({moltype:[molid]})
		# Attribute the molecule to an existing group
		else:
			if moltype in groups.keys():
				groups.get(moltype).append(molid)
			else:
				groups.update({moltype:[molid]})
	return groups
# Function fitplane should be used to give the equation of a plane that determined precisely by three points or "close" to more than three points
# "close" here means that each cosine of angle between normal vector and the vector from the center of cluster of points to one point should be as close to 90 degrees (perpendicular) as possible
# Function fitplane return four parameters: a, b, c, d. The equation of wanted plane is "ax+by+cz+d=0"
def fitplane(pts):
	# Use leastsq module from scipy to converge the least squared summation
    from scipy.optimize import leastsq
	# Function planeq is used for leastsq
	# p should be a normal vector(to initialize) or objects can be converted to a vector
    def planeq(p, x, y, z):
        if isinstance(p, np.ndarray) and p.size == 3 and p.ndim == 1:
            pass
        elif isinstance(p, (tuple, list)) and len(p) == 3 and all([isinstance(v, (int, float)) for v in p]):
            p = np.array(p)
        else:
            raise TypeError("Unsupported parameter given for plane equation.")
        pt = np.array([x, y, z])
		# CAN NOT use np.dot(pt, pt) to replace (x**2 + y**2 + z**2) - They are NOT equal
        return (np.dot(p, pt))/np.sqrt((p.dot(p) * (x**2 + y**2 + z**2)))
    if isinstance(pts, (tuple, list)) and len(pts) > 2 and all([isinstance(pt, np.ndarray) and pt.size == 3 and pt.ndim == 1 for pt in pts]):
		# Get the normal vector from the first three points
        normal = np.cross(pts[1]-pts[0], pts[2]-pts[0])
        a, b, c = normal
        if len(pts) == 3:
			# When only three points are given, return the precise equation of the plane
            d = -(normal*pts[0]).sum()
            return a, b, c, d
        else:
			# When more than three points are given, just return equation of a "close" plane
			# adjpts is a N*3 matrix that each row equals to a vector from the center of N points to the point
            adjpts = np.array(pts) - np.average(np.array(pts), axis=0)
			# Transpose of adjpts is a 3*N matrix, three rows cooresponds to x, y, z, respectively
            x, y, z = adjpts.transpose()
            a, b, c = leastsq(planeq, normal, args=(x, y, z))[0]
            d = -np.array([a, b, c]).dot(np.average(np.array(pts), axis=0))
            return a, b, c, d
    else:
        raise ValueError("Unsupported parameter given for plane function in function fitplane from module gromacs.")
def rotate(vector, axis, cosang, sinang):
	x, y, z = axis
	rotmatrix = (1-cosang)*np.outer(axis, axis) + sinang*np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]]) + cosang*np.eye(3)
	return np.dot(rotmatrix, vector)
# Function hbond only judge whether a hydrogen bond exists with these three coordinates
def hbond(donor, hydro, acceptor, distance = 0.35, angle = 30):
	if not isinstance(distance, float) and not isinstance(angle, (int, float)) and distance > 0 and 0 <= angle <= 180:
		raise ValueError("Wrong parameters. Function hbond expects distance(nm) > 0 and 0 <= angle(degree) <= 180.")
	dist = math.sqrt(np.dot(donor-acceptor, donor-acceptor))
	ang = math.degrees(math.acos(np.dot(hydro-donor, acceptor-donor)/math.sqrt(np.dot(hydro-donor, hydro-donor)*np.dot(acceptor-donor, acceptor-donor))))
	if dist < distance and ang < angle:
		return True
	else:
		return False
# Parameters molA and molB both should be pack of mol id, mol name, atom names and mol coordinate matrix
def hbondmol(molA, molB):
	def hbondlist(moltype):
		# Attention: atom ids in each list should start from 0, not 1
		# Add new molecular type here when try to analysis a new system
		if moltype == "BCD":
			donlist = [(0, 1), (4, 5), (10, 11), (19, 20), (22, 23), (25, 26), (33, 34), (36, 37), (39, 40), (47, 48), (50, 51), (53, 54), (61, 62), (64, 65), (67, 68), (75, 76), (78, 79), (81, 82), (87, 88), (90, 91), (96, 97)]
			acclist = [0, 4, 7, 10, 12, 14, 17, 19, 22, 25, 28, 31, 33, 36, 39, 42, 45, 47, 50, 53, 56, 59, 61, 64, 67, 70, 73, 75, 78, 81, 84, 87, 90, 93, 96]
		elif moltype == "C8A":
			donlist = [(2, 11)]
			acclist = [0, 2]
		elif moltype == "C8O":
			donlist = [(8, 9)]
			acclist = [8]
		elif moltype == "C8N":
			donlist = [(8, 9), (8, 10)]
			acclist = [8]
		elif moltype == "QAC":
			donlist = []
			acclist = []
		elif moltype == "NO3":
			donlist = []
			acclist = [1, 2, 3]
		elif moltype == "SOL":
			donlist = [(0, 1), (0, 2)]
			acclist = [0]
		else:
			raise ValueError("Molecular type %s has no corresponding data term. Please check the original code and append data record to function hbondlist." % moltype)
		return donlist, acclist
	# Unpack molA and molB
	mola_id, mola_name, mola_atom, mola_coord = molA
	molb_id, molb_name, molb_atom, molb_coord = molB
	# Lists of donor and acceptor atom ids
	donlista, acclista = hbondlist(mola_name)
	donlistb, acclistb = hbondlist(molb_name)
	# Initialize log
	log = []
	# Count total hbond number between molA and molB
	hbondn = 0
	for (a, h) in donlista:
		for b in acclistb:
			# Detail(mol id + mol name + atom name)
			don = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name + " " + '{0:>4}'.format(mola_atom[a]))
			hyd = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name + " " + '{0:>4}'.format(mola_atom[h]))
			acc = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name + " " + '{0:>4}'.format(molb_atom[b]))
			# Coordinate
			donor = mola_coord[a]
			hydro = mola_coord[h]
			acceptor = molb_coord[b]
			if hbond(donor, hydro, acceptor) is True:
				hbondn += 1
				dist = math.sqrt(np.dot(donor-acceptor, donor-acceptor))
				ang = math.degrees(math.acos(np.dot(hydro-donor, acceptor-donor)/math.sqrt(np.dot(hydro-donor, hydro-donor)*np.dot(acceptor-donor, acceptor-donor))))
				log.append(don + " | " + hyd + " | " + acc + " | " + '{0:>9.3f}'.format(dist) + " | " + '{0:>6.1f}'.format(ang))
	for b in acclista:
		for (a,h) in donlistb:
			# Detail(mol id + mol name + atom name)
			don = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name + " " + '{0:>4}'.format(molb_atom[a]))
			hyd = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name + " " + '{0:>4}'.format(molb_atom[h]))
			acc = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name + " " + '{0:>4}'.format(mola_atom[b]))
			# Coordinate
			donor = molb_coord[a]
			hydro = molb_coord[h]
			acceptor = mola_coord[b]
			if hbond(donor, hydro, acceptor) is True:
				hbondn += 1
				dist = math.sqrt(np.dot(donor-acceptor, donor-acceptor))
				ang = math.degrees(math.acos(np.dot(hydro-donor, acceptor-donor)/math.sqrt(np.dot(hydro-donor, hydro-donor)*np.dot(acceptor-donor, acceptor-donor))))
				log.append(don + " | " + hyd + " | " + acc + " | " + '{0:>9.3f}'.format(dist) + " | " + '{0:>6.1f}'.format(ang))
	return hbondn, log
class Moledit(object):
	def __init__(self, pdbfile):
		self.mols = pdb2molmatrix(pdbfile)
		self.parentdir = os.path.split(pdbfile)[0]
		self.pdbname = os.path.split(pdbfile)[1]
		self.command()
	def command(self):
		print("Molecules Editing for GROMACS")
		cmd = input(">>> ").strip()
		lcmd = {"move":self.move, "moveto":self.moveto, "rotate":self.orrient, "random-rotate":self.random_orrient,"write":self.write}
		# Read commands from console
		while cmd != "exit":
			if cmd.split(" ")[0] not in lcmd.keys():
				print("[Info] Unknown command.")
			else:
				try:
					lcmd.get(cmd.split(" ")[0])(*([convert_input_type(x) for x in cmd.split(" ")[1:]]))
				except Exception as e:
					print("[Error] "+str(e)+".")
			cmd = input(">>> ").strip()
		save_and_exit = input("Save modifications?(Y/N)").strip().upper()
		while save_and_exit not in ("Y", "N"):
			save_and_exit = input("Save modifications?(Y/N)").strip().upper()
		if save_and_exit == "Y":
			self.write()
		input("Press any key to exit...")
	def cut(self, centermol = 1, shape = "sphere", parameter = (1.0,)):
		# Check
		if centermol not in self.mols.keys():
			raise ValueError("Can not find target molecule in the file provided.")
		# Calculate coordinate of target molecular center
		center = np.average(self.mols.get(centermol).get("coordinate").transpose(), axis=1)
		# Cut a sphere around target molecular center
		if shape in ["sphere"]:
			inside, outside = [], []
			radius = parameter[0]
			# Scan all molecules
			for i in range(1, len(self.mols)+1):
				mol = self.mols.get(i)
				molcenter = np.average(mol.get("coordinate").transpose(), axis=1)
				# Put a molecule "inside" when its center is within the sphere
				if np.dot(molcenter-np.array(center), molcenter-np.array(center)) < radius**2:
					inside.update({i:mol})
				# Or put it "outside"
				else:
					outside.append({i:mol})
			return inside, outside
		else:
			raise ValueError("Unsupported operation name given.")
	def move(self, movemol = 1, vector = [0., 0., 0.]):
		# Check
		if isinstance(vector, np.ndarray) and vector.ndim == 1 and vector.size == 3:
			pass
		elif isinstance(vector, (tuple, list)) and len(vector) == 3 and all([isinstance(v, (int, float)) for v in vector]):
			vector = np.array(vector)
		else:
			raise TypeError("Get wrong parameters. \"vector\" should be able to convert to a 1*3 matrix.")
		if movemol not in self.mols.keys():
			raise ValueError("Can not find target molecule in the file provided.")
		# Move
		self.mols.get(movemol).update({"coordinate":self.mols.get(movemol).get("coordinate")+vector})
	def moveto(self, movemol = 1, point = [0., 0., 0.]):
		# Check
		if isinstance(point, np.ndarray) and point.ndim == 1 and point.size == 3:
			pass
		elif isinstance(point, (tuple, list)) and len(point) == 3 and all([isinstance(v, (int, float)) for v in point]):
			point = np.array(point)
		else:
			raise TypeError("Get wrong parameters. \"point\" should be able to convert to a 1*3 matrix.")
		if movemol not in self.mols.keys():
			raise ValueError("Can not find target molecule in the file provided.")
		# Move target molecule to destination
		# The molecule goes to where its center overlap the point
		center = np.average(self.mols.get(movemol).get("coordinate").transpose(), axis=1)
		self.mols.get(movemol).update({"coordinate":self.mols.get(movemol).get("coordinate")+point-center})
	def orrient(self, movemol = 1, vector = [0., 0., 1.]):
		# Check
		if isinstance(vector, np.ndarray) and vector.ndim == 1 and vector.size == 3:
			pass
		elif isinstance(vector, (tuple, list)) and len(vector) == 3 and all([isinstance(v, (int, float)) for v in vector]):
			vector = np.array(vector)
		else:
			raise TypeError("Get wrong parameters. \"vector\" should be able to convert to a 1*3 matrix.")
		if all([v == 0 for v in vector]):
			raise ValueError("Get wrong parameters. \"vector\" should not be a zero vector.")
		if movemol not in self.mols.keys():
			raise ValueError("Can not find target molecule in the file provided.")
		# Normalize
		vector = vector/math.sqrt(vector.dot(vector))
		# Coordinate matrix of target molecule
		molcoord = self.mols.get(movemol).get("coordinate")
		center = np.average(molcoord, axis=0)
		# Get the normal vector determined by molecular coordinates
		plps = fitplane([v for v in molcoord])
		normal = np.array(plps[0:3])
		normal = normal/math.sqrt(normal.dot(normal))
		# Rotation axis
		axis = np.cross(normal, vector)
		axis = axis/math.sqrt(axis.dot(axis))
		# Rotation angle
		cosang = np.dot(normal, vector)
		sinang = math.sqrt(1-cosang**2)
		# Rotate
		newcoord = []
		for i in range(len(molcoord)):
			newcoord.append(rotate(molcoord[i]-center, axis, cosang, sinang) + center)
		self.mols.get(movemol).update({"coordinate":np.array(newcoord)})
	def random_orrient(self, movemol):
		import random
		self.orrient(movemol, np.array([random.random() for i in range(3)]))
	def write(self, copy = False):
		if copy is True:
			pdbname = input("Please enter the name to save as:")
		else:
			pdbname = self.pdbname
		try:
			molmatrix2pdb(self.mols, os.path.join(self.parentdir, pdbname))
		except Exception as e:
			print("[Error] "+str(e)+".")
		else:
			print("[Info] File %s has been written successfully." % pdbname)
class TrajAnalysis(object):
	def load(self, dirpath, num, name, frametype = "pdb"):
		# Initialize
		self.trajprof = {}
		self.trajcoord = []
		self.trajn = num
		# Load trajectory profile from the first frame(should end with number 0)
		firstframe = pdb2molmatrix(os.path.join(dirpath, name + "0" + "." + frametype))
		# Remove coordinate and velocity terms
		# Add atom id term for each atom, begin counting
		# Atom ids would be count from 1
		atomn = 1
		for molid in firstframe.keys():
			mol = firstframe.get(molid)
			mol_atomn = len(mol.get("atom"))
			molprof = {"name":mol.get("name"), "atom":mol.get("atom"), "atomid":[atomn+i for i in range(mol_atomn)]}
			atomn += mol_atomn
			self.trajprof.update({molid:molprof})
		# atomn equals to total counting of atoms plus initial 1
		self.atomn = atomn - 1
		# Only load coordinates data of each frame (including the first one) into a matrix (N*3)
		for i in range(num):
			mols = pdb2molmatrix(os.path.join(dirpath, name + repr(i) + "." + frametype))
			self.trajcoord.append(np.vstack([mol.get("coordinate") for mol in mols.values()]))
			print("[Info] File %s is loaded successfully." % (name + repr(i) + "." + frametype))
	def distance(self, atomA, atomB):
		# Check
		if not isinstance(atomA, int) or not ((isinstance(atomB, int) or (isinstance(atomB, (tuple, list)) and all([isinstance(v, int) for v in atomB])))):
			raise ValueError("Wrong parameters. Function distance expects the first parameter should be ID for central atom and the second one be either one ID or some IDs packed in a tuple or list for surrounding atom(s).")
		if not (isinstance(boxpara, (tuple, list)) and len(boxpara) == 3 and all([isinstance(v, (int, float)) for v in boxpara])) or not isinstance(pbc, bool):
			raise ValueError("Wrong parameters. Function distance expects default parameter boxpara be a tuple os list with three number each be parameter of box in one dimension.")
		if atomA > self.atomn or not (((isinstance(atomB, int) and atomB <= self.atomn) or (isinstance(atomB, (tuple, list)) and all([x <= self.atomn for x in atomB])))):
			raise ValueError("Can not find pair atom(s) in the matrix provided.")
		# Format
		if isinstance(atomB, (tuple, list)):
			titleline = '{0:>6}'.format("frames") + " | " + " | ".join(['{0:>5}'.format(atomA) + "-" + '{0:>5}'.format(atomBi) for atomBi in atomB])
		else:
			titleline = '{0:>6}'.format("frames") + " | " + '{0:>5}'.format(atomA) + "-" + '{0:>5}'.format(atomB)
		# Initialize log and counting
		log = []
		framen = -1
		for framecoord in self.trajcoord:
			framen += 1
			coordA = framecoord[atomA-1]
			if isinstance(atomB, (tuple, list)):
				dists = []
				for atomBi in atomB:
					# Original Coordinate
					pair0 = coordA - framecoord[atomBi-1]
					pairs = [pair0]
					findist = min([math.sqrt(x.dot(x)) for x in pairs])
					dists.append(findist)
				log.append('{0:>6}'.format(framen) + " | " + " | ".join(['{0:>11.3f}'.format(d) for d in dists]))
			else:
				# Original Coordinate
				pair0 = coordA - framecoord[atomBi-1]
				pairs = [pair0]
				findist = min([math.sqrt(x.dot(x)) for x in pairs])
				log.append('{0:>6}'.format(framen) + " | " + '{0:>11.3f}'.format(findist))
		# Print to screen
		print(titleline)
		for line in log:
			print(line)
	def dist2plane(self, atom, plane):
		# Check
		if not isinstance(atom, int) or not (isinstance(plane, (tuple, list)) and len(plane) >= 3 and all([isinstance(v, int) for v in plane])):
			raise ValueError("Wrong parameters. Function distance expects two parameters: the former should be ID for central atom and the latter should be atom IDs in the reference plane.")
		if not (isinstance(boxpara, (tuple, list)) and len(boxpara) == 3 and all([isinstance(v, (int, float)) for v in boxpara])) or not isinstance(pbc, bool):
			raise ValueError("Wrong parameters. Function distance expects default parameter boxpara be a tuple os list with three number each be parameter of box in one dimension.")
		if atom > self.atomn or all([x <= self.atomn for x in plane]) is False:
			raise ValueError("Can not find atom(s) in the matrix provided.")
		# Format
		titleline = '{0:>6}'.format("frames") + " | " + '{0:>8}'.format("distance")
		# Initialize log and counting
		log = []
		framen = 0
		for framecoord in self.trajcoord:
			framen += 1
			# Original Coordinate
			coord = framecoord[atom-1]
			coords= [coord]
			# Function fitplane only accepts a package of at least three points' coordinate (better near a plane) in the form of either tuple or list, each element of the package should be a 1-dim 3-size numpy.ndarray object
			plps = fitplane([framecoord[x-1] for x in plane])	# Short of plane parameters
			# Attention: plps should be a tuple with 4 elements, denoted as a, b, c, d. The plane equation should be ax+by+cz=0
			plnormal = np.array(plps[0:3])
			pldist = min([abs(plnormal.dot(c) + plps[3])/(plnormal.dot(plnormal))**(1/2) for c in coords])
			log.append('{0:>6}'.format(framen) + " | " + '{0:>8.3f}'.format(pldist))
		# Print to screen
		print(titleline)
		for line in log:
			print(line)
	def hbondgrp(self, moltypeA, moltypeB):
		def hbondmolgrp(prof, framecoord, moltypeA, moltypeB):
			# Grouping
			grp = cluster(prof)
			# Initialize log and counting
			log = []
			hbondn = []
			molids = []
			for mola_id in grp.get(moltypeA):
				for molb_id in grp.get(moltypeB):
					mola = prof.get(mola_id)
					molA = (
						mola_id,
						mola.get("name"),
						mola.get("atom"),
						np.vstack([framecoord[l-1] for l in mola.get("atomid")])
					)
					molb = prof.get(molb_id)
					molB = (
						molb_id,
						molb.get("name"),
						molb.get("atom"),
						np.vstack([framecoord[l-1] for l in molb.get("atomid")])
					)
					mol_hbondn, mol_log = hbondmol(molA, molB)
					hbondn.append(mol_hbondn)
					log.extend(mol_log)
					if mol_hbondn > 0:
						molids.extend([mola_id, molb_id])
			molids = list(set(molids))
			molids.sort()
			return hbondn, log, molids
		# Check
		profgrp = cluster(self.trajprof)
		if moltypeA not in profgrp.keys() or moltypeB not in profgrp.keys():
			raise ValueError("Wrong parameters. Function hbondgrp expects at least two parameters, both of them should be three-letter tags of corresponding molecules.")
		# Format
		frametitle = '{0:>14}'.format("donor") + " | " + '{0:>14}'.format("hydro") + " | " + '{0:>14}'.format("acceptor") + " | " + '{0:>9}'.format("distance") + " | " + '{0:>6}'.format("angle")
		# Start counting number of frames
		framen = 0
		# Initialize dict
		tothbondn = {}
		hbonds = {}
		for framecoord in self.trajcoord:
			hbondn, log, molids = hbondmolgrp(self.trajprof, framecoord, moltypeA, moltypeB)
			frame_hbondn = sum(hbondn)
			if frame_hbondn > 0:
				print("<frame %d>" % framen)
				print(frametitle)
				for line in log:
					print(line)
			tothbondn.update({framen:frame_hbondn})
			hbonds.update({framen:molids})
			# Counting of frames start from 0
			framen += 1
		print('{0:>6}'.format("frame") + " | " + '{0:>6}'.format("hbondn"))
		for i in range(self.trajn):
			print('{0:>6}'.format(i) + " | " + '{0:>6}'.format(tothbondn.get(i)))
		return tothbondn, hbonds
	# Extract some molecules out of one frame
	def write(self, frameid, molid):
		molid = list(set(molid))
		molid.sort()
		if isinstance(frameid, int) and frameid < self.trajn and ((isinstance(molid, int) and molid in self.trajprof.keys()) or (isinstance(molid, (tuple, list)) and all([v in self.trajprof.keys() for v in molid]))):
			molmatrix = {}
			framecoord = self.trajcoord[frameid]
			moln = 0
			# Reconstruct the molecular matrix
			for i in molid:
				moln += 1
				molprof = self.trajprof.get(i)
				mol = {"name":molprof.get("name"), "atom":molprof.get("atom")}
				mol.update({"coordinate":np.vstack([framecoord[l-1] for l in molprof.get("atomid")])})
				molmatrix.update({moln:mol})
			molmatrix2pdb(molmatrix).write("newframe%d.pdb" %frameid)
		else:
			raise ValueError("Wrong parameters. Function write expects two parameters, the former should be id for the frame and the latter should be either an id for the molecular matrix or a tuple or list for several ids.")
