import os
import math
import random
import numpy as np

from gmx.inputFunctions import convert_input_type
from gmx.io.log import *
from gmx.structure.molMatrix import *
from gmx.math.common import cluster
from gmx.math.common import fitplane
from gmx.math.common import rotate
from gmx.math.hbond import hbondmol

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
		moln = 0
		for molid in firstframe.keys():
			mol = firstframe.get(molid)
			mol_atomn = len(mol.get("atom"))
			molprof = {"name":mol.get("name"), "atom":mol.get("atom"), "atomid":[atomn+i for i in range(mol_atomn)]}
			atomn += mol_atomn
			moln += 1
			self.trajprof.update({molid:molprof})
		# atomn equals to total counting of atoms plus initial 1
		self.atomn = atomn - 1
		# moln equals to total molecule number
		self.moln = moln
		# Only load coordinates data of each frame (including the first one) into a matrix (N*3)
		for i in range(num):
			mols = pdb2molmatrix(os.path.join(dirpath, name + repr(i) + "." + frametype))
			self.trajcoord.append(np.vstack([mol.get("coordinate") for mol in mols.values()]))
			print("[Info] File %s is loaded successfully." % (name + repr(i) + "." + frametype))
	def distance(self, atomA, atomB):
		# Check
		if not isinstance(atomA, int) or not ((isinstance(atomB, int) or (isinstance(atomB, (tuple, list)) and all([isinstance(v, int) for v in atomB])))):
			raise ValueError("Wrong parameters. Function distance expects the first parameter should be ID for central atom and the second one be either one ID or some IDs packed in a tuple or list for surrounding atom(s).")
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
	def moldist(self, molid, atomid):
		# Check
		if not isinstance(molid, int) or not ((isinstance(atomid, int) or (isinstance(atomid, (tuple, list)) and all([isinstance(v, int) for v in atomid])))):
			raise ValueError("Wrong parameters. Function distance expects the first parameter should be ID for the molecule and the second one be either one ID or some IDs packed in a tuple or list for surrounding atom(s).")
		if molid > self.moln or not (((isinstance(atomid, int) and atomid <= self.atomn) or (isinstance(atomid, (tuple, list)) and all([x <= self.atomn for x in atomid])))):
			raise ValueError("Can not find pair atom(s) in the matrix provided.")
		# Format
		if isinstance(atomid, (tuple, list)):
			titleline = '{0:>6}'.format("frames") + " | " + " | ".join(['{0:>5}'.format(molid) + "-" + '{0:>5}'.format(atomi) for atomi in atomid])
		else:
			titleline = '{0:>6}'.format("frames") + " | " + '{0:>5}'.format(molid) + "-" + '{0:>5}'.format(atomid)
		# Initialize log and counting
		log = []
		framen = -1
		for framecoord in self.trajcoord:
			framen += 1
			# Calculate molecular center coordinate
			molatomids = self.trajprof.get(molid).get("atomid")
			center = np.average(framecoord[min(molatomids)-1:max(molatomids)], axis=0)
			if isinstance(atomid, (tuple, list)):
				dists = []
				for atomi in atomid:
					# Original Coordinate
					pair0 = center - framecoord[atomi-1]
					pairs = [pair0]
					findist = min([math.sqrt(x.dot(x)) for x in pairs])
					dists.append(findist)
				log.append('{0:>6}'.format(framen) + " | " + " | ".join(['{0:>11.3f}'.format(d) for d in dists]))
			else:
				# Original Coordinate
				pair0 = center - framecoord[atomid-1]
				pairs = [pair0]
				findist = min([math.sqrt(x.dot(x)) for x in pairs])
				log.append('{0:>6}'.format(framen) + " | " + '{0:>11.3f}'.format(findist))
		# Print to screen
		print(titleline)
		for line in log:
			print(line)
	def atomLocation(self, atomid):
		# Check
		if (
		isinstance(atomid, int) and atomid <= self.atomn
		):
			# Format
			titleline = '{0:>8}'.format("Frame") + " | " + '{0:>6}'.format("X") + " | " + '{0:>6}'.format("Y") + " | " + '{0:>6}'.format("Z")
			# Initialize log and counting
			log = []
			framen = -1
			sum = np.array([0., 0., 0.])
			for framecoord in self.trajcoord:
				framen += 1
				# Calculate molecular center coordinate
				location = framecoord[atomid-1]
				sum += location
				log.append('{0:>8}'.format(framen) + " | " + " | ".join(['{0:>6.3f}'.format(c) for c in location]))
			# Print to screen
			print(titleline)
			for line in log:
				print(line)
			print('{0:>8}'.format("Average") + " | " + " | ".join(['{0:>6.3f}'.format(c / (framen + 1)) for c in sum]))
		elif (
		isinstance(atomid, (tuple, list)) and all([isinstance(a, int) and a <= self.atomn for a in atomid])
		):
			# Format
			titleline = '{0:>8}'.format("Frame") + "|" + " | ".join(['{0:>5}'.format(a) + "X" + " | " + '{0:>5}'.format(a) + "Y" + " | " + '{0:>5}'.format(a) + "Z" for a in atomid])
			# Initialize log and counting
			log = []
			framen = -1
			sum = np.array([0.] * (3 * len(atomid)))
			for framecoord in self.trajcoord:
				framen += 1
				location = []
				# Calculate molecular center coordinate
				for a in atomid:
					location.extend(framecoord[a-1].tolist())
				location = np.array(location)
				sum += location
				log.append('{0:>8}'.format(framen) + " | " + " | ".join(['{0:>6.3f}'.format(c) for c in location]))
			# Print to screen
			print(titleline)
			for line in log:
				print(line)
			print('{0:>8}'.format("Average") + " | " + " | ".join(['{0:>6.3f}'.format(c / (framen + 1)) for c in sum]))
		else:
			raise ValueError("Wrong parameters. Function atomLocation expects an integer parameter which is the selected atom ID, or a tuple or list each of which is atom ID wanted.")
	def molcenter(self, molid):
		# Check
		if not isinstance(molid, int):
			raise ValueError("Wrong parameters. Function molcenter expects an integer parameter which is the selected molecule ID.")
		if molid > self.moln:
			raise ValueError("Can not find pair molecule in the matrix provided.")
		# Format
		titleline = '{0:>8}'.format("Frame") + " | " + '{0:>6}'.format("X") + " | " + '{0:>6}'.format("Y") + " | " + '{0:>6}'.format("Z")
		# Initialize log and counting
		log = []
		framen = -1
		sum = np.array([0., 0., 0.])
		for framecoord in self.trajcoord:
			framen += 1
			# Calculate molecular center coordinate
			molatomids = self.trajprof.get(molid).get("atomid")
			center = np.average(framecoord[min(molatomids)-1:max(molatomids)], axis=0)
			sum += center
			log.append('{0:>8}'.format(framen) + " | " + " | ".join(['{0:>6.3f}'.format(c) for c in center]))
		# Print to screen
		print(titleline)
		for line in log:
			print(line)
		print('{0:>8}'.format("Average") + " | " + " | ".join(['{0:>6.3f}'.format(c / (framen + 1)) for c in sum]))
	def dist2plane(self, atom, plane):
		# Check
		if not isinstance(atom, int) or not (isinstance(plane, (tuple, list)) and len(plane) >= 3 and all([isinstance(v, int) for v in plane])):
			raise ValueError("Wrong parameters. Function distance expects two parameters: the former should be ID for central atom and the latter should be atom IDs in the reference plane.")
		if atom > self.atomn or all([x <= self.atomn for x in plane]) is False:
			raise ValueError("Can not find atom(s) in the matrix provided.")
		# Format
		titleline = '{0:>6}'.format("frames") + " | " + '{0:>8}'.format("distance")
		# Initialize log and counting
		log = []
		framen = -1
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
	def moldist2plane(self, molid, plane):
		# Check
		if not isinstance(molid, int) or not (isinstance(plane, (tuple, list)) and len(plane) >= 3 and all([isinstance(v, int) for v in plane])):
			raise ValueError("Wrong parameters. Function distance expects two parameters: the former should be ID for the molecule and the latter should be atom IDs in the reference plane.")
		if molid > self.moln or all([x <= self.atomn for x in plane]) is False:
			raise ValueError("Can not find atom(s) in the matrix provided.")
		# Format
		titleline = '{0:>6}'.format("frames") + " | " + '{0:>8}'.format("distance")
		# Initialize log and counting
		log = []
		framen = -1
		for framecoord in self.trajcoord:
			framen += 1
			# Calculate molecular center coordinate
			molatomids = self.trajprof.get(molid).get("atomid")
			center = np.average(framecoord[min(molatomids)-1:max(molatomids)], axis=0)
			coords= [center]
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
		print('{0:>6}'.format("frames") + " | " + '{0:>6}'.format("hbondn"))
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
