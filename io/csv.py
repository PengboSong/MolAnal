from gmx.other.input_func import convert_input_type, is_int, is_float

def readline_csv(line, split_symbol = ';', comment_symbol = '#'):
	if not line.startswith(comment_symbol):
		return line.split(comment_symbol, 1)[0].strip().split(split_symbol)
	else:
		return None

def readline_csv_spec(line, row):
	assert isinstance(row, int), "Invalid row number."
	assert isinstance(line, str), "Can not parse the given line at row %d." % row
	# Newline; Oldline; AtomName; ResType; ResId; Occup; TempFactor; SegName
    split_line = readline_csv(line)
	if isinstance(split_line, list) and len(split_line) == 8:
		new_id, old_id, atom_name, res_type, res_id, occup, temp_factor, segment = split_line
		assert is_int(new_id), "Term Newline at row %d column 1 is not an integer." % row
		new_id = int(new_id)
		if new_id != 0:
			assert is_int(old_id), "Term Oldline at row %d column 2 is not an integer." % row
			old_id = int(old_id)

			atom_name = convert_input_type(atom_name)
			assert isinstance(atom_name, str), "Term AtomName at row %d column 3 is not characters." % row

			res_type = convert_input_type(res_type)
			assert isinstance(res_type, str), "Term ResType at row %d column 4 is not characters." % row

			assert is_int(res_id), "Term ResId at row %d column 5 is not an integer." % row
			res_id = int(res_id)

			assert is_float(occup), "Term Occup at row %d column 6 is not a real number." % row
			occup = float(occup)

			assert is_float(temp_factor), "Term TempFactor at row %d column 7 is not a real number." % row
			occup = float(occup)

			segment = convert_input_type(segment)
			assert isinstance(segment, str), "Term SegName at row %d column 8 is not characters." % row

			return {
				"new_atom_serial_number":new_id,
				"atom_serial_number":old_id,
				"atom_name":atom_name,
				"residual_name":res_type,
				"residue_sequence_number":res_id,
				"occupancy":occup,
				"temperature_factor":temp_factor,
				"element_symbol":segment
			}
		else:
			return None
	else:
		return None
