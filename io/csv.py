from gmx.inputFunctions import convert_input_type

def readline_csv(line, split_symbol = ';', comment_symbol = '#'):
	if not line.startswith(comment_symbol):
		return line.split(comment_symbol, 1)[0].strip().split(split_symbol)
	else:
		return None

def readline_csv_spec(line, row):
	# Newline; Oldline; AtomName; ResType; ResId; Occup; TempFactor; SegName
    splitLine = readline_csv(line)
	if isinstance(splitLine, list) and len(splitLine) == 8:
		new_id, old_id, atom_name, res_type, res_id, occup, temp_factor, segment = splitLine
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
