import os, sys
import math
import random
import numpy as np
import time

from gmx.inputFunctions import convert_input_type
from gmx.io.log import *
from gmx.structure.molMatrix import *
from gmx.math.common import cluster
from gmx.math.common import fitplane
from gmx.math.common import rotate
from gmx.math.common import convertSeconds
from gmx.math.hbond import hbondmol
from gmx.math.rmsd import rmsdMol

from custom.general import setProgressBar

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
	def __init__(self):
		self.trajprof = {}
		self.trajcoord = []
		self.trajn = 0
		self.atomn = 0
		self.moln = 0

		self._savePath = "."
		self._loadPath = "."

	def getSavePath(self):
		return self._savePath
	def setSavePath(self, path):
		if os.path.isdir(path):
			self._savePath = path
		else:
			print("[Error] Invalid path. The given path should be a existing directory.")

	def getLoadPath(self):
		return self._loadPath
	def setLoadPath(self, path):
		if os.path.isdir(path):
			self._loadPath = path
		else:
			print("[Error] Invalid path. The given path should be a existing directory.")

	def clearAll(self):
		self.trajprof = {}
		self.trajcoord = []
		self.trajn = 0
		self.trajstartn = 0
		self.trajendn = 0
		self.atomn = 0
		self.moln = 0

	def listCommand(self):
		cmdList = [x for x in dir(self) if not x.startswith("_") and callable(getattr(self, x))]
		cmdList.sort()
		for i in range(len(cmdList)):
			print(str(i + 1) + ": " + cmdList[i])

	def command(self):
		def exit():
			exitConfirm = input("Exit?(Y/N)").strip().upper()
			while exitConfirm not in ("Y", "N"):
				exitConfirm = input("Save modifications?(Y/N)").strip().upper()
			if exitConfirm == "Y":
				return True
			else:
				return False

		print("Trajectory Analysis for GROMACS")
		cmd = input(">>> ").strip()
		# Read commands from console
		while True:
			inputCmdList = cmd.split()
			if inputCmdList[0] == "exit":
				if exit():
					break
			elif not hasattr(self, inputCmdList[0]):
				print("[Info] Invalid command.")
			else:
				try:
					getattr(self, inputCmdList[0])(*([convert_input_type(x) for x in inputCmdList[1:]]))
				except Exception as e:
					print("[Error] "+str(e)+".")
			cmd = input(">>> ").strip()
		input("Press any key to exit...")

	def load(self, start, end, name, frametype = "pdb"):
		def framePath(index):
			if isinstance(index, int):
				frameName = name + str(index) + "." + frametype
				return os.path.join(self._loadPath, frameName)
			else:
				return None

		# Initialize
		self.trajprof = {}
		self.trajcoord = []
		self.trajn = end - start + 1
		self.trajstartn = start
		self.trajendn = end
		# Load trajectory profile from the first frame(should end with number 0)
		firstframe = pdb2molmatrix(framePath(start))
		# Remove coordinate and velocity terms
		# Add atom id term for each atom, begin counting
		# Atom ids would be count from 1
		atomn = 1
		moln = 0
		for molid in firstframe.keys():
			mol = firstframe.get(molid)
			mol_atomn = len(mol.get("atom"))
			molprof = {
			"name":mol.get("name"),
			"atom":mol.get("atom"),
			"atomid":[atomn+i for i in range(mol_atomn)]
			}
			atomn += mol_atomn
			moln += 1
			self.trajprof.update({molid:molprof})
		# atomn equals to total counting of atoms minus initial 1
		self.atomn = atomn - 1
		# moln equals to total molecule number
		self.moln = moln
		# Only load coordinates data of each frame (including the first one) into a matrix (N*3)
		sys.stdout.write("[Info] Start reading files.\n")
		for i in range(self.trajstartn, self.trajendn + 1):
			mols = pdb2molmatrix(framePath(i))
			self.trajcoord.append(np.vstack([mol.get("coordinate") for mol in mols.values()]))
			sys.stdout.write(setProgressBar(100, round((i / self.trajn) * 100)))
			sys.stdout.flush()
		sys.stdout.write(setProgressBar(100, 100))
		sys.stdout.write("\n[Info] All File have loaded successfully.\r\n")

	def writeLog(self, name, log, suffix = ".txt"):
		if not isinstance(name, str) or not (isinstance(log, (tuple, list)) and all([isinstance(l, str) for l in log])):
			raise TypeError("[Error] Wrong parameters. Function writeLog expects the first parameter to be a string and the second one to be a string list (or tuple).")
		else:
			# yyyymmdd-HHMMSS
			timeFormat = time.strftime("%Y%m%d_%H%M%S")
			fileName = name + "_" + timeFormat + suffix
			with open(os.path.join(self._savePath, fileName), "w") as f:
				f.writelines(line + '\n' for line in log)

	def checkAtomID(self, atomID):
		if isinstance(atomID, int) and atomID <= self.atomn:
			return True
		else:
			return False
	def checkAtomIDList(self, atomIDList):
		if isinstance(atomIDList, (tuple, list)) and all([self.checkAtomID(a) for a in atomIDList]):
			return True
		else:
			return False

	def checkMolID(self, molID):
		if isinstance(molID, int) and molID <= self.moln:
			return True
		else:
			return False
	def checkMolIDList(self, molIDList):
		if isinstance(molIDList, (tuple, list)) and all([self.checkMolID(m) for m in molIDList]):
			return True
		else:
			return False
	def checkFrameID(self, frameID):
		if isinstance(frameID, int) and self.trajstartn <= frameID <= self.trajendn:
			return True
		else:
			return False
	def checkFrameIDList(self, frameIDList):
		if isinstance(frameIDList, (tuple, list)) and all([self.checkMolID(f) for f in frameIDList]):
			return True
		else:
			return False

	def molCenter(self, frameID, molID):
		atomIDList = self.trajprof.get(molID).get("atomid")
		framecoord = self.trajprof[frameID - self.trajstartn]
		center = np.average(framecoord[min(atomIDList) - 1 : max(atomIDList)], axis=0)
		return center
	def molCenterWithCoord(self, framecoord, molID):
		atomIDList = self.trajprof.get(molID).get("atomid")
		center = np.average(framecoord[min(atomIDList) - 1 : max(atomIDList)], axis=0)
		return center

	def distance(self, atomA, atomB):
		# Check
		if self.checkAtomID(atomB):
			atomB = tuple(atomB)
		if not self.checkAtomID(atomA) or not self.checkAtomIDList(atomB):
			raise ValueError("Wrong parameters. Function distance expects the first parameter to be the central atom ID and the second one to be an ID list (or tuple) for pair atom(s).")

		# Title
		titleline = '{0:>6}'.format("frames") + " | " + " | ".join(['{0:>5}'.format(atomA) + "-" + '{0:>5}'.format(b) for b in atomB])
		log = [titleline]
		# Start Counting
		# Be aware that counting number plus 1 at the beginning of the loop
		# So it should begin from the former integer of start frame ID
		# The same below
		framen = self.trajstartn - 1
		for framecoord in self.trajcoord:
			framen += 1
			coordA = framecoord[atomA-1]
			dists = []
			for b in atomB:
				pair0 = coordA - framecoord[b - 1]
				findist = math.sqrt(pair0.dot(pair0))
				dists.append(findist)
			log.append('{0:>6}'.format(framen) + " | " + " | ".join(['{0:>11.3f}'.format(d) for d in dists]))

		self.writeLog("distance", log)

	def distanceMol(self, molid, atomid):
		# Check
		if self.checkAtomID(atomid):
			atomid = tuple(atomid)
		if not self.checkMolID(molid) or not self.checkAtomIDList(atomid):
			raise ValueError("Wrong parameters. Function distance expects the first parameter to be the selected molecular ID and the second one to be an ID list (or tuple) for pair atom(s).")

		# Title
		titleline = '{0:>6}'.format("frames") + " | " + " | ".join(['{0:>5}'.format(molid) + "-" + '{0:>5}'.format(a) for a in atomid])
		log = [titleline]
		# Start Counting
		framen = self.trajstartn - 1
		for framecoord in self.trajcoord:
			framen += 1
			center = self.molCenterWithCoord(framecoord, molid)
			dists = []
			for a in atomid:
				pair0 = center - framecoord[a - 1]
				findist = math.sqrt(pair0.dot(pair0))
				dists.append(findist)
			log.append('{0:>6}'.format(framen) + " | " + " | ".join(['{0:>11.3f}'.format(d) for d in dists]))

		self.writeLog("distance_mol", log)

	def location(self, atomid):
		# Check
		if self.checkAtomID(atomid):
			atomid = tuple(atomid)
		if not self.checkAtomIDList(atomid):
			raise ValueError("Wrong parameters. Function atomLocation expects one parameter to be an ID list (or tuple) for selected atom(s).")

		# Title
		titleline = '{0:>8}'.format("Frame") + " | " + " | ".join(['{0:>5}'.format(a) + "X" + " | " + '{0:>5}'.format(a) + "Y" + " | " + '{0:>5}'.format(a) + "Z" for a in atomid])
		log = [titleline]
		# Start Counting
		framen = self.trajstartn - 1
		sum = np.zeros((3 * len(atomid),))
		for framecoord in self.trajcoord:
			framen += 1
			locations = []
			for a in atomid:
				locations.extend(framecoord[a-1].tolist())
			locations = np.array(locations)
			sum += locations
			log.append('{0:>8}'.format(framen) + " | " + " | ".join(['{0:>6.3f}'.format(c) for c in locations]))
		log.append('{0:>8}'.format("Average") + " | " + " | ".join(['{0:>6.3f}'.format(c / (framen + 1)) for c in sum]))

		self.writeLog("location", log)

	def locationMol(self, molid):
		# Check
		if self.checkMolID(molid):
			molid = tuple(molid)
		if not self.checkMolIDList(molid):
			raise ValueError("Wrong parameters. Function molLocation expects one parameter to be an ID list (or tuple) for selected molecule(s).")

		# Title
		titleline = '{0:>8}'.format("Frame") + " | " + " | ".join(['{0:>5}'.format(m) + "X" + " | " + '{0:>5}'.format(m) + "Y" + " | " + '{0:>5}'.format(m) + "Z" for m in molid])
		log = [titleline]
		# Start Counting
		framen = self.trajstartn - 1
		sum = np.zeros((3 * len(molid),))
		for framecoord in self.trajcoord:
			framen += 1
			centers = []
			for m in molid:
				centers.extend(self.molCenterWithCoord(framecoord, m))
			centers = np.array(centers)
			sum += centers
			log.append('{0:>8}'.format(framen) + " | " + " | ".join(['{0:>6.3f}'.format(c) for c in centers]))
		log.append('{0:>8}'.format("Average") + " | " + " | ".join(['{0:>6.3f}'.format(c / (framen + 1)) for c in sum]))

		self.writeLog("location_mol", log)

	def dist2plane(self, atomid, plane):
		# Check
		if not self.checkAtomID(atomid) or not self.checkAtomIDList(plane):
			raise ValueError("Wrong parameters. Function dist2plane expects the first parameter to be the central atom ID and the second one to be an ID list (or tuple) for atom(s) in the reference plane.")
		else:
			if len(plane) < 3:
				raise ValueError("Need at least 3 points to fit a plane.")

		# Format
		titleline = '{0:>6}'.format("frames") + " | " + '{0:>8}'.format("distance")
		log = [titleline]
		# Start Counting
		framen = self.trajstartn - 1
		for framecoord in self.trajcoord:
			framen += 1
			coord = framecoord[atomid - 1]
			# Function fitplane only accepts a package of at least three points' coordinate (better near a plane) in the form of either tuple or list, each element of the package should be a 1-dim 3-size numpy.ndarray object
			plps = fitplane([framecoord[x - 1] for x in plane])	# Short of plane parameters
			# Attention: plps should be a tuple with 4 elements, denoted as a, b, c, d. The plane equation should be ax+by+cz=0
			plnormal = np.array(plps[0:3])
			pldist = abs(plnormal.dot(coord) + plps[3])/math.sqrt(plnormal.dot(plnormal))
			log.append('{0:>6}'.format(framen) + " | " + '{0:>8.3f}'.format(pldist))

		self.writeLog("dist2plane", log)

	def dist2planeMol(self, molid, plane):
		# Check
		if not self.checkMolID(molid) or not self.checkAtomIDList(plane):
			raise ValueError("Wrong parameters. Function dist2planeMol expects the first parameter to be the selected molecular ID and the second one to be an ID list (or tuple) for atom(s) in the reference plane.")
		else:
			if len(plane) < 3:
				raise ValueError("Need at least 3 points to fit a plane.")

		# Title
		titleline = '{0:>6}'.format("frames") + " | " + '{0:>8}'.format("distance")
		log = [titleline]
		# Start Counting
		framen = self.trajstartn - 1
		for framecoord in self.trajcoord:
			framen += 1
			center = self.molCenterWithCoord(framecoord, molid)
			# Function fitplane only accepts a package of at least three points' coordinate (better near a plane) in the form of either tuple or list, each element of the package should be a 1-dim 3-size numpy.ndarray object
			plps = fitplane([framecoord[x - 1] for x in plane])	# Short of plane parameters
			# Attention: plps should be a tuple with 4 elements, denoted as a, b, c, d. The plane equation should be ax+by+cz=0
			plnormal = np.array(plps[0:3])
			pldist = abs(plnormal.dot(center) + plps[3])/math.sqrt(plnormal.dot(plnormal))
			log.append('{0:>6}'.format(framen) + " | " + '{0:>8.3f}'.format(pldist))

		self.writeLog("dist2plane_mol", log)

	def rmsd(self, framen):
		# Check
		if not self.checkFrameID(framen):
			raise ValueError("Wrong parameters. Function rmsd expects one parameter to be the selected frame ID.")

		# Title
		titleLine = '{0:>6}'.format("frames") + " | " + '{0:>8}'.format("RMSD")
		log = [titleLine]
		# Reference frame
		refMatrix = self.trajcoord[framen - self.trajstartn]
		weightFactor = np.ones(self.atomn)
		# Start Counting
		framen = self.trajstartn
		for i in range(self.trajn):
			framen += 1
			thisMatrix = self.trajcoord[i]
			rmsdv = rmsdMol(refMatrix, thisMatrix, weightFactor)
			log.append('{0:>6}'.format(framen) + " | " + '{0:>8.3f}'.format(rmsdv))

		self.writeLog("rmsd_mol", log)

	def clusterSingleLinkage(self, cutoff):
		# Construct RMSD Matrix
		rmsdMatrix = np.zeros((self.trajn, self.trajn))
		weightFactor = np.ones(self.atomn)
		startTime = time.time()
		sys.stdout.write("[Info] Start calculating RMSD matrix.\n")
		for i in range(self.trajn - 1):
			for j in range(i + 1, self.trajn):
				rmsdMatrix[i][j] = rmsdMol(self.trajcoord[i], self.trajcoord[j], weightFactor)
				rmsdMatrix[j][i] = rmsdMatrix[i][j]
			sys.stdout.write(setProgressBar(100, round((i / self.trajn) * 100)))
			sys.stdout.flush()
		endTime = time.time()
		sys.stdout.write(setProgressBar(100, 100))
		sys.stdout.write("\n[Info] RMSD matrix calculation normally ends.\r\n")
		print("[Info] Total time usage = %s." % convertSeconds(endTime - startTime))

		rmsdLog = [";".join(['{0:>6.3f}'.format(i) for i in r]) for r in rmsdMatrix]

		# Cluster
		clusters = [[0 + self.trajstartn]]
		for i in range(1, self.trajn):
			newGroupFlag = True
			for group in clusters:
				# rmsdMatrix index equals to frame index minus start frame index, start from 0
				if all([rmsdMatrix[i][g - self.trajstartn] < cutoff for g in group]):
					group.append(i + self.trajstartn)
					newGroupFlag = False
					break
			if newGroupFlag is True:
				clusters.append([i + self.trajstartn])

		titleLine = '{0:>6}'.format("Groupn") + " | " + '{0:>24}'.format("Group members (Framen)")
		log = [titleLine]
		groupn = 0
		for group in clusters:
			groupn += 1
			log.append('{0:>6}'.format(groupn) + " | " + ', '.join([str(g) for g in group]))

		self.writeLog("cluster_single_linkage", log)
		self.writeLog("cluster_rmsd_matrix", rmsdLog, ".csv")

		return rmsdMatrix

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
		# Title
		frametitle = '{0:>14}'.format("donor") + " | " + '{0:>14}'.format("hydro") + " | " + '{0:>14}'.format("acceptor") + " | " + '{0:>9}'.format("distance") + " | " + '{0:>6}'.format("angle")
		log = []
		# Start Counting
		# Counting at the end of the loop
		framen = self.trajstartn
		# Initialize dict
		tothbondn = {}
		hbonds = {}
		for framecoord in self.trajcoord:
			hbondn, orglog, molids = hbondmolgrp(self.trajprof, framecoord, moltypeA, moltypeB)
			frame_hbondn = sum(hbondn)
			if frame_hbondn > 0:
				log.append("<frame %d>" % framen)
				log.append(frametitle)
				log.extend(orglog)
			tothbondn.update({framen:frame_hbondn})
			hbonds.update({framen:molids})
			framen += 1
		log.append('{0:>6}'.format("frames") + " | " + '{0:>6}'.format("hbondn"))
		for i in range(self.trajstartn, self.trajendn + 1):
			log.append('{0:>6}'.format(i) + " | " + '{0:>6}'.format(tothbondn.get(i)))
		return tothbondn, hbonds

	# Extract some molecules out of one frame
	def writeNewFrame(self, frameid, molid):
		molid = list(set(molid))
		molid.sort()
		if not self.checkFrameID(frameid) or not self.checkMolIDList(molid):
			raise ValueError("Wrong parameters. Function writeNewFrame expects two parameters, the former one to be the selected frame ID and the latter one to be an ID list (or tuple) for the wanted molecule(s).")
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
		molmatrix2pdb(molmatrix, "newframe%d.pdb" % frameid)

	def writeGroups(self, pairs, filename):
		# pairs should be a tuple list
		# each tuple has two elements, the former one should be frame ID, and the latter one should be molecular ID
		if (isinstance(pairs, list) and all([isinstance(p, tuple) and len(p) == 2 and self.checkFrameID(p[0]) and self.checkMolIDList(p[1]) for p in pairs])):
			molmatrix = {}
			moln = 0
			for pair in pairs:
				frameid = pair[0]
				molid = pair[1]
				for i in molid:
					moln += 1
					molprof = self.trajprof.get(i)
					mol = {"name":molprof.get("name"), "atom":molprof.get("atom")}
					mol.update({"coordinate":np.vstack([self.trajcoord[frameid][l-1] for l in molprof.get("atomid")])})
					molmatrix.update({moln:mol})
			molmatrix2pdb(molmatrix, "%s.pdb" % filename)
		else:
			raise ValueError("Wrong parameters. Function writeGroups expects two parameters, the former one to be a tuple list in which each element should be a pack of frame ID and molecular ID list and the latter one to be the wanted file name without suffix.")
