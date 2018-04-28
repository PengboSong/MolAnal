from gmx.inputFunctions import convert_input_type

def readline_gro(line, row):
    rowText = format(row, 'd')
	# Read by column number
	# Molecular id: 1-5, character, right
	mol_id = convert_input_type(line[0:5].strip())
	if not isinstance(mol_id, int):
		raise ValueError("Wrong data type. Molecular id at row " + rowText + " column 1-5 is not an integer.")
	# Molecular name: 6-8, character, left
	mol_name = convert_input_type(line[5:8].strip())
	if not isinstance(mol_name, str):
		raise ValueError("Wrong data type. Molecular name at row " + rowText + " column 6-8 is not characters.")
	# Chain identifier: 10, character
	chain = convert_input_type(line[9:10].strip())
	if chain and not isinstance(chain, str):
		raise ValueError("Wrong data type. Chain identifier at row " + rowText + " column 10 is not a character.")
	# Atom name: 12-15, character, right
	atom_name = convert_input_type(line[11:15].strip())
	if not isinstance(atom_name, str):
		raise ValueError("Wrong data type. Atom name at row " + rowText + " column 12-15 is not characters.")
	# Atom id: 16-20, integer, right
	atom_id = convert_input_type(line[15:20].strip())
	if not isinstance(atom_id, int):
		raise ValueError("Wrong data type. Atom id at row " + rowText + " column 16-20 is not an integer.")
	# Coordinate X(A): 21-28, real(8:3), right
	x = convert_input_type(line[20:28].strip())
	if not isinstance(x, float):
		raise ValueError("Wrong data type. Coordinate X at row " + rowText + " column 21-28 is not a real number.")
	# Coordinate Y(A): 29-36, real(8:3), right
	y = convert_input_type(line[28:36].strip())
	if not isinstance(y, float):
		raise ValueError("Wrong data type. Coordinate Y at row " + rowText + " column 29-36 is not a real number.")
	# Coordinate Z(A): 37-44, real(8:3), right
	z = convert_input_type(line[36:44].strip())
	if not isinstance(z, float):
		raise ValueError("Wrong data type. Coordinate Z at row " + rowText + " column 37-44 is not a real number.")
	if len(line.strip()) > 44:
		# Velocity Vx: 45-52, read(8:4), right
		vx = convert_input_type(line[44:52].strip())
		if not isinstance(vx, float):
			raise ValueError("Wrong data type. Velocity Vx at row " + rowText + " column 45-52 is not a real number.")
		# Velocity Vy: 53-60, read(8:4), right
		vy = convert_input_type(line[52:60].strip())
		if not isinstance(vy, float):
			raise ValueError("Wrong data type. Velocity Vy at row " + rowText + " column 53-60 is not a real number.")
		# Velocity Vz: 61-68, read(8:4), right
		vz = convert_input_type(line[60:68].strip())
		if not isinstance(vz, float):
			raise ValueError("Wrong data type. Velocity Vz at row " + rowText + " column 61-68 is not a real number.")
	else:
		vx, vy, vz = 0., 0., 0.
	return {"mol_id":mol_id, "mol_name":mol_name, "chain":chain, "atom_name":atom_name, "atom_id":atom_id, "x":x, "y":y, "z":z, "vx":vx, "vy":vy, "vz":vz}
