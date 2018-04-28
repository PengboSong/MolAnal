import math
import numpy as np

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
	elif moltype == "DCF":
		donlist = [(8, 22)]
		acclist = [18, 19]
	elif moltype == "LAS":
		donlist = []
		acclist = [7, 8, 9]
	elif moltype == "SOL":
		donlist = [(0, 1), (0, 2)]
		acclist = [0]
	else:
		raise ValueError("Molecular type %s has no corresponding data term. Please check the original code and append data record to function hbondlist." % moltype)
	return donlist, acclist

# Parameters molA and molB both should be pack of mol id, mol name, atom names and mol coordinate matrix
def hbondmol(molA, molB):
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
