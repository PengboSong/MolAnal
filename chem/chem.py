def num2element(num):
    elements = ["H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "CL", "KR"]
    if num > len(elements):
        return None
    else:
        return elements[num - 1]

def element2num(element):
	elements = {"H":1, "HE":2, "LI":3, "BE":4, "B":5, "C":6, "N":7, "O":8, "F":9, "NE":10, "NA":11, "MG":12, "AL":13, "SI":14, "P":15, "S":16, "CL":17, "AR":18, "K":19, "CA":20, "SC":21, "TI":22, "V":23, "CR":24, "MN":25, "FE":26, "CO":27, "NI":28, "CU":29, "ZN":30, "GA":31, "GE":32, "AS":33, "SE":34, "CL":35, "KR":36}
	if element in elements.keys():
		return elements.get(element)
	else:
		return None

def element2mass(element):
	elements = {"H":1.00794, "HE":4.0026, "LI":6.941, "BE":9.01218, "B":10.811, "C":12.0107, "N":14.0067, "O":15.9994, "F":18.9984, "NE":20.1797, "NA":22.9898, "MG":24.3050, "AL":26.9815, "SI":28.0855, "P":30.9738, "S":32.065, "CL":35.453, "AR":39.948, "K":39.0983, "CA":40.078, "SC":44.9559, "TI":47.867, "V":50.9415, "CR":51.9961, "MN":54.9380, "FE":55.845, "CO":58.9332, "NI":58.6934, "CU":63.546, "ZN":65.409, "GA":69.723, "GE":72.64, "AS":74.9216, "SE":78.96, "CL":79.904, "KR":39.948}
	if element in elements.keys():
		return elements.get(element)
	else:
		return None
