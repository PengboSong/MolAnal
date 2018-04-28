from gmx.inputFunctions import convert_input_type

def readline_pdb(line, row):
    rowText = format(row, 'd')
	# Read by column number
	# Record type: 1-6, character, left
	record = convert_input_type(line[0:6].strip())
	if not isinstance(record, str):
		raise ValueError("Wrong data type. Record type at row " + rowText + " column 1-6 is not characters.")
	# Atom serial number: 7-11, integer, right
	atom_id = convert_input_type(line[6:11].strip())
	if not isinstance(atom_id, int):
		raise ValueError("Wrong data type. Atom serial number at row " + rowText + " column 7-11 is not an integer.")
	# Atom name: 13-16, character, left
	atom_name = convert_input_type(line[12:16].strip())
	if not isinstance(atom_name, str):
		raise ValueError("Wrong data type. Atom name at row " + rowText + " column 13-16 is not characters.")
	# Alternate location indicator: 17, character
	location = convert_input_type(line[16:17].strip())
	if location and not isinstance(location, str):
		raise ValueError("Wrong data type. Alternate location at row " + rowText + " column 17 is not a character.")
	# Residual name: 18-20, character, right
	res_name = convert_input_type(line[17:20].strip())
	if res_name and not isinstance(res_name, str):
		raise ValueError("Wrong data type. Residual name at row " + rowText + " column 18-20 is not characters.")
	# Chain identifier: 22, character
	chain = convert_input_type(line[21:22].strip())
	if chain and not isinstance(chain, str):
		raise ValueError("Wrong data type. Chain identifier at row " + rowText + " column 22 is not characters.")
	# Residue sequence number: 23-26, interger, right
	res_id = convert_input_type(line[22:26].strip())
	if res_id and not isinstance(res_id, int):
		raise ValueError("Wrong data type. Residue sequence number at row " + rowText + " column 23-26 is not an integer.")
	# Code for insertions of residuals: 27, character
	code = convert_input_type(line[26:27].strip())
	if code and not isinstance(code, str):
		raise ValueError("Wrong data type. Code for insertions of residuals at row " + rowText + " column 27 is not a character.")
	# Coordinate X(A): 31-38, real(8:3), right
	x = convert_input_type(line[30:38].strip())
	if not isinstance(x, float):
		raise ValueError("Wrong data type. Coordinate X at row " + rowText + " column 31-38 is not a real number.")
	# Coordinate Y(A): 39-46, real(8:3), right
	y = convert_input_type(line[38:46].strip())
	if not isinstance(y, float):
		raise ValueError("Wrong data type. Coordinate Y at row " + rowText + " column 39-46 is not a real number.")
	# Coordinate Z(A): 47-54, real(8:3), right
	z = convert_input_type(line[46:54].strip())
	if not isinstance(z, float):
		raise ValueError("Wrong data type. Coordinate Z at row " + rowText + " column 47-54 is not a real number.")
	# Occupancy: 55-60, real(6:2), right
	occup = convert_input_type(line[54:60].strip())
	if occup and not isinstance(occup, float):
		raise ValueError("Wrong data type. Occupancy at row " + rowText + " column 55-60 is not a real number.")
	# Temperature factor: 61-66, real(6:2), right
	temp_factor = convert_input_type(line[60:66].strip())
	if temp_factor and not isinstance(temp_factor, float):
		raise ValueError("Wrong data type. Temperature factor at row " + rowText + " column 61-66 is not a real number.")
	# Segment identifier: 73-76, character, left
	segment = convert_input_type(line[72:76].strip())
	if segment and not isinstance(segment, str):
		raise ValueError("Wrong data type. Segment identifier at row " + rowText + " column 73-76 is not characters.")
	# Element symbol: 77-78, character, right
	element = convert_input_type(line[76:78].strip())
	if element and not isinstance(element, str):
		raise ValueError("Wrong data type. Element symbol at row " + rowText + " column 77-78 is not characters.")
	return {"record_type":record, "atom_serial_number":atom_id, "atom_name":atom_name, "alternate_location_identifier":location, "residual_name":res_name, "chain_identifier":chain, "residue_sequence_number":res_id, "code_for_insertions_of_residues":code, "coordinate_x":x, "coordinate_y":y, "coordinate_z":z, "occupancy":occup, "temperature_factor":temp_factor, "segment_identifier":segment, "element_symbol":element}
