import collections

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

def is_int(s):
	assert isinstance(s, str)
	is_int_flag = False
	try:
		int(s)
	except Exception:
		pass
	else:
		is_int_flag = True
	return is_int_flag

def is_float(s):
	assert isinstance(s, str)
	is_float_flag = False
	if is_int(s) is False:
		try:
			float(s)
		except Exception:
			pass
		else:
			is_float_flag = True
	return is_float_flag

def split_skip_domain(s, split_symbol = ','):
	assert isinstance(s, str)
	left_paired_symbol = {'(':')', '[':']', '{':'}'}

	paired_symbol = ''
	paired_list = [0]
	left_domain_count = 0
	right_domain_count = 0
	find_paired_symbol = True
	for i, char in enumerate(s):
		if find_paired_symbol is False:
			if char == paired_symbol and left_domain_count == right_domain_count:
				paired_list.append(i + 1)
				left_domain_count = 0
				right_domain_count = 0
				find_paired_symbol = True
			elif char in left_paired_symbol.keys():
				left_domain_count += 1
			elif char in left_paired_symbol.values():
				right_domain_count += 1
		elif char in left_paired_symbol.keys() and find_paired_symbol:
			paired_symbol = left_paired_symbol.get(s[i])
			paired_list.append(i)
			tmp_domain_id = i
			find_paired_symbol = False
	paired_list.append(len(s))

	assert find_paired_symbol is True, "Find incomplete domain."
	result = []
	for i in range(len(paired_list)):
		if i % 2 == 0:
			split_part = s[paired_list[i]:paired_list[i + 1]].split(split_symbol)
			if result:
				result[-1] += split_part[0]
				result.extend(split_part[1:])
			else:
				result.extend(split_part)
		elif i + 1 < len(paired_list):
			protect_part = s[paired_list[i]:paired_list[i + 1]]
			if result:
				result[-1] += protect_part
			else:
				result.append(protect_part)
	return result

def convert_input_type(string):
	result = None
	# Only handle strings
	if isinstance(string, str):
		string = string.strip()
		# First, judge whether the string is an integer
		if is_int(string):
			result = int(string)
		# Second, judge whether it is a floating number
		elif is_float(string):
			result = float(string)
		elif string.startswith('[') and string.endswith(']'):
			content = split_skip_domain(string[1:-1], ',')
			result = [convert_input_type(c) for c in content]
		elif string.startswith('(') and string.endswith(')'):
			content = split_skip_domain(string[1:-1], ',')
			result = tuple([convert_input_type(c) for c in content])
		elif string.startswith('{') and string.endswith('}'):
			content = split_skip_domain(string[1:-1], ',')
			if all([':' in c for c in content]) is True:
				result = {}
				for c in content:
					key, value = c.split(':', 1)
					key = convert_input_type(key)
					value = convert_input_type(value)
					if isinstance(key, collections.Hashable):
						result.update({key: value})
			elif any([':' in c for c in content]) is False:
				result = set([convert_input_type(c) for c in content])
			else:
				result = string
		else:
			result = string
	return result
