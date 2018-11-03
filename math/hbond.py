import math
import numpy as np

from gmx.chem.molecules import hbondlist, hbondlist_direction
from gmx.math.common import cluster, fitplane, rotate

hbond_distance = 0.35
hbond_degree = 30
hbond_dir_degree = 30

def hbond(donor, hydro, acceptor):
	dist = math.sqrt(np.dot(donor-acceptor, donor-acceptor))
	ang = math.degrees(math.acos(np.dot(hydro-donor, acceptor-donor)/math.sqrt(np.dot(hydro-donor, hydro-donor)*np.dot(acceptor-donor, acceptor-donor))))
	return dist, ang

def hbond_dir_sp2_withnormal(normal, center, neighbor):
	v = neighbor - center
	unitv = v / math.sqrt(v.dot(v))
	rv1 = rotate(unitv, normal, 1/2, math.sqrt(3)/2)
	rv2 = rotate(unitv, normal, -1/2, math.sqrt(3)/2)
	return rv1, rv2

def hbond_dir_sp2_withoutnormal(center, neighbor1, neighbor2):
	v1 = neighbor1 - center
	v2 = neighbor2 - center
	v = v1 + v2
	unitv = v / math.sqrt(v.dot(v))
	return unitv

def hbond_dir_sp3_twoadjacent(center, neighbor1, neighbor2):
	v1 = neighbor1 - center
	v2 = neighbor2 - center
	arbaxis = v1 + v2
	axis = arbaxis / math.sqrt(arbaxis.dot(arbaxis))
	rv1 = rotate(unitv, axis, 0, 1)
	rv2 = rotate(unitv, axis, 0, 1)
	return -rv1, -rv2

def hbond_dir_sp3_threeadjacent(center, neighbor1, neighbor2, neighbor3):
	v1 = neighbor1 - center
	v2 = neighbor2 - center
	v3 = neighbor3 - center
	v = v1 + v2 + v3
	unitv = v / math.sqrt(v.dot(v))
	return unitv

def hbond_dir(coord, dirlistAll):
	result = []
	for dirlist in dirlistAll:
		spn = dirlist[0]
		nbn = len(dirlist) - 1
		if spn == "sp2":
			if nbn == 2:
				normal = np.array(fitplane([coord[k] for k in range(coord.shape[0])]))[0:3]
				v1, v2 = hbond_dir_sp2_withnormal(normal, *[coord[i] for i in dirlist[1:]])
				result.append(v1)
				result.append(v2)
			elif nbn == 3:
				v = hbond_dir_sp2_withoutnormal(*[coord[i] for i in dirlist[1:]])
				result.append(v)
		elif spn == "sp3":
			if nbn == 3:
				v1, v2 = hbond_dir_sp3_twoadjacent(*[coord[i] for i in dirlist[1:]])
				result.append(v1)
				result.append(v2)
			elif nbn == 4:
				v = hbond_dir_sp3_threeadjacent(*[coord[i] for i in dirlist[1:]])
				result.append(v)
	if len(result) > 0:
		return result
	else:
		return None

def verifyHbond(
	mola_id, mola_name, mola_atom, mola_coord, donlista, acclista,
	molb_id, molb_name, molb_atom, molb_coord, donlistb, acclistb,
	hbondn, log
):
	for (a, h) in donlista:
		for b in acclistb:
			# Coordinate
			donor = mola_coord[a]
			hydro = mola_coord[h]
			acceptor = molb_coord[b]
			dist, ang = hbond(donor, hydro, acceptor)
			if dist < hbond_distance and ang < hbond_degree:
				hbondn += 1
				# Detail(mol id + mol name + atom name)
				don = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name) + " " + '{0:>4}'.format(mola_atom[a])
				hyd = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name) + " " + '{0:>4}'.format(mola_atom[h])
				acc = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name) + " " + '{0:>4}'.format(molb_atom[b])
				log.append(don + " | " + hyd + " | " + acc + " | " + '{0:>9.3f}'.format(dist) + " | " + '{0:>6.1f}'.format(ang))
	for b in acclista:
		for (a, h) in donlistb:
			# Coordinate
			donor = molb_coord[a]
			hydro = molb_coord[h]
			acceptor = mola_coord[b]
			dist, ang = hbond(donor, hydro, acceptor)
			if dist < hbond_distance and ang < hbond_degree:
				hbondn += 1
				# Detail(mol id + mol name + atom name)
				don = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name) + " " + '{0:>4}'.format(molb_atom[a])
				hyd = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name) + " " + '{0:>4}'.format(molb_atom[h])
				acc = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name) + " " + '{0:>4}'.format(mola_atom[b])
				log.append(don + " | " + hyd + " | " + acc + " | " + '{0:>9.3f}'.format(dist) + " | " + '{0:>6.1f}'.format(ang))
	return hbondn, log

def verifyHbond_dir(
	mola_id, mola_name, mola_atom, mola_coord, donlista, acclista, dirlista,
	molb_id, molb_name, molb_atom, molb_coord, donlistb, acclistb, dirlistb,
	hbondn, log, directiondb
):
	for (a, h) in donlista:
		for b in acclistb:
			# Coordinate
			donor = mola_coord[a]
			hydro = mola_coord[h]
			acceptor = molb_coord[b]
			dist, ang = hbond(donor, hydro, acceptor)
			if dist < hbond_distance and ang < hbond_degree:
				withDirection = True

				# Check whether directions of molB(acceptor) are already calculated
				if molb_id in directiondb.keys():
					directionb = directiondb.get(molb_id)
				else:
					# If not, calculate them and add these values to the dict
					directionb = hbond_dir(molb_coord, dirlistb)
					if directionb:
						directiondb.update({molb_id: directionb})
					else:
						withDirection = False

				isHydroBond = False
				ang_ds = []
				# Consider whether the hydrogen bond is close to the most possible directions
				if withDirection is True:
					for d in directionb:
						ang_d = math.degrees(math.acos(np.dot(hydro-acceptor, d)/math.sqrt(np.dot(hydro-acceptor, hydro-acceptor))))
						if ang_d < hbond_dir_degree:
							isHydroBond = True
							ang_ds.append(ang_d)
							break
				# When without valid values, only consider distance(usually <0.35 nm) and angle(usually <30 degrees)
				else:
					isHydroBond = True

				if len(ang_ds) == 0:
					ang_d = 0
				else:
					ang_d = min(ang_ds)

				if isHydroBond is True:
					hbondn += 1
					# Detail(mol id + mol name + atom name)
					don = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name) + " " + '{0:>4}'.format(mola_atom[a])
					hyd = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name) + " " + '{0:>4}'.format(mola_atom[h])
					acc = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name) + " " + '{0:>4}'.format(molb_atom[b])
					log.append(don + " | " + hyd + " | " + acc + " | " + '{0:>9.3f}'.format(dist) + " | " + '{0:>6.1f}'.format(ang) + ", " + '{0:>6.1f}'.format(ang_d))

	for b in acclista:
		for (a, h) in donlistb:
			# Coordinate
			donor = molb_coord[a]
			hydro = molb_coord[h]
			acceptor = mola_coord[b]
			dist, ang = hbond(donor, hydro, acceptor)
			if dist < hbond_distance and ang < hbond_degree:
				withDirection = True

				# Check whether directions of molA(acceptor) are already calculated
				if mola_id in directiondb.keys():
					directiona = directiondb.get(mola_id)
				else:
					# If not, calculate them and add these values to the dict
					directiona = hbond_dir(mola_coord, dirlista)
					if directiona:
						directiondb.update({mola_id: directiona})
					else:
						withDirection = False

				isHydroBond = False
				ang_ds = []
				# Consider whether the hydrogen bond is close to the most possible directions
				if withDirection is True:
					for d in directiona:
						ang_d = math.degrees(math.acos(np.dot(hydro-acceptor, d)/math.sqrt(np.dot(hydro-acceptor, hydro-acceptor))))
						if ang_d < hbond_dir_degree:
							isHydroBond = True
							ang_ds.append(ang_d)
							break
				# When without valid values, only consider distance(usually <0.35 nm) and angle(usually <30 degrees)
				else:
					isHydroBond = True

				if len(ang_ds) == 0:
					ang_d = 0
				else:
					ang_d = min(ang_ds)

				if isHydroBond is True:
					hbondn += 1
					# Detail(mol id + mol name + atom name)
					don = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name) + " " + '{0:>4}'.format(molb_atom[a])
					hyd = '{0:>5}'.format(molb_id) + " " + '{0:>4}'.format(molb_name) + " " + '{0:>4}'.format(molb_atom[h])
					acc = '{0:>5}'.format(mola_id) + " " + '{0:>4}'.format(mola_name) + " " + '{0:>4}'.format(mola_atom[b])
					log.append(don + " | " + hyd + " | " + acc + " | " + '{0:>9.3f}'.format(dist) + " | " + '{0:>6.1f}'.format(ang) + ", " + '{0:>6.1f}'.format(ang_d))

	return hbondn, log, directiondb

# Parameters molA and molB both should be pack of mol id, mol name, atom names and mol coordinate matrix
def hbondmol(molA, molB, log):
	# Unpack molA and molB
	mola_id, mola_name, mola_atom, mola_coord = molA
	molb_id, molb_name, molb_atom, molb_coord = molB
	# Lists of donor and acceptor atom ids
	donlista, acclista = hbondlist(mola_name)
	donlistb, acclistb = hbondlist(molb_name)

	# Initialize counting times
	hbondn = 0

	return verifyHbond(
		mola_id, mola_name, mola_atom, mola_coord, donlista, acclista,
		molb_id, molb_name, molb_atom, molb_coord, donlistb, acclistb,
		hbondn, log
	)

def hbondmol_multiple(molA, molb_name, molb_num, molb_ids, molb_atoms, molb_coords, log):
	# molA info
	mola_id, mola_name, mola_atom, mola_coord = molA
	donlista, acclista = hbondlist(mola_name)

	# molB info
	donlistb, acclistb = hbondlist(molb_name)

	hbondns = []
	for ni in range(molb_num):
		molb_id = molb_ids[ni]
		molb_atom = molb_atoms[ni]
		molb_coord = molb_coords[ni]

		# Initialize counting times for each molB
		hbondn = 0
		hbondn, log = verifyHbond(
			mola_id, mola_name, mola_atom, mola_coord, donlista, acclista,
			molb_id, molb_name, molb_atom, molb_coord, donlistb, acclistb,
			hbondn, log
		)
		hbondns.append(hbondn)
	return hbondns, log

def hbondmol_dir(molA, molB, log, directiondb):
	# Unpack molA and molB
	mola_id, mola_name, mola_atom, mola_coord = molA
	molb_id, molb_name, molb_atom, molb_coord = molB
	# Lists of donor and acceptor atom ids
	donlista, acclista = hbondlist(mola_name)
	dirlista = hbondlist_direction(mola_name)
	donlistb, acclistb = hbondlist(molb_name)
	dirlistb = hbondlist_direction(molb_name)
	# Initialize counting times
	hbondn = 0
	return verifyHbond_dir(
		mola_id, mola_name, mola_atom, mola_coord, donlista, acclista, dirlista,
		molb_id, molb_name, molb_atom, molb_coord, donlistb, acclistb, dirlistb,
		hbondn, log, directiondb
	)

def hbondmol_dir_multiple(molA, molb_name, molb_num, molb_ids, molb_atoms, molb_coords, log, directiondb):
	# molA info
	mola_id, mola_name, mola_atom, mola_coord = molA
	donlista, acclista = hbondlist(mola_name)
	dirlista = hbondlist_direction(mola_name)

	# molB info
	donlistb, acclistb = hbondlist(molb_name)
	dirlistb = hbondlist_direction(molb_name)

	hbondns = []
	for ni in range(molb_num):
		molb_id = molb_ids[ni]
		molb_atom = molb_atoms[ni]
		molb_coord = molb_coords[ni]

		# Initialize counting times for each molB
		hbondn = 0
		hbondn, log, directiondb = verifyHbond_dir(
			mola_id, mola_name, mola_atom, mola_coord, donlista, acclista, dirlista,
			molb_id, molb_name, molb_atom, molb_coord, donlistb, acclistb, dirlistb,
			hbondn, log, directiondb
		)
		hbondns.append(hbondn)
	return hbondns, log, directiondb

def hbondmolgrp(prof, framecoord, moltypeA, moltypeB, lpdir = False):
	# Grouping
	grp = cluster(prof)

	# Initialize log and countings
	log = []
	hbondns = []
	molids = []
	if lpdir is True:
		directiondb = {}

	grpA = grp.get(moltypeA)
	grpB = grp.get(moltypeB)
	molan = len(grpA)
	molbn = len(grpB)

	if molan > molbn:
		molan, molbn = molbn, molan
		grpA, grpB = grpB, grpA

	for mola_id in grpA:
		mola = prof.get(mola_id)
		molA = (
			mola_id,
			mola.get("name"),
			mola.get("atom"),
			np.vstack([framecoord[l-1] for l in mola.get("atomid")])
		)
		if 10 * molan < molbn:
			oneMolb = prof.get(grpB[0])
			molb_name = oneMolb.get("name")
			molb_num = len(grpB)
			molb_ids = grpB
			molb_atoms = [prof.get(id).get("atom") for id in grpB]
			molb_coords = [np.vstack([framecoord[l-1] for l in prof.get(id).get("atomid")]) for id in grpB]
			if lpdir is True:
				hbondns_single, log, directiondb = hbondmol_dir_multiple(molA, molb_name, molb_num, molb_ids, molb_atoms, molb_coords, log, directiondb)
			else:
				hbondns_single, log = hbondmol_multiple(molA, molb_name, molb_num, molb_ids, molb_atoms, molb_coords, log)
			hbondns.extend(hbondns_single)
			hbondnsArray = np.array(hbondns_single)
			grpbArray = np.array(grpB)
			if np.any(hbondnsArray > 0) is True:
				molids.append(mola_id)
			molids.extend(grpbArray[hbondnsArray > 0].tolist())
		else:
			for molb_id in grpB:
				molb = prof.get(molb_id)
				molB = (
					molb_id,
					molb.get("name"),
					molb.get("atom"),
					np.vstack([framecoord[l-1] for l in molb.get("atomid")])
				)
				if lpdir is True:
					hbondn, log, directiondb = hbondmol_dir(molA, molB, log, directiondb)
				else:
					hbondn, log = hbondmol(molA, molB, log)
				hbondns.append(hbondn)
				if hbondn > 0:
					molids.extend([mola_id, molb_id])
	molids = list(set(molids))
	molids.sort()
	return hbondns, log, molids
