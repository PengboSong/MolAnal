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
