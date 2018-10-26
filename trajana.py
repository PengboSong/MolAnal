import os, time, math
import numpy as np

from gmx.structure.molMatrix import molmatrix2pdb, pdb2molmatrix
from gmx.structure.rmsdMatrix import checkRmsdMatrix, upperRmsdMatrix, reshapeRmsdMatrix

from gmx.math.clust import newCluster, linkStructure, renumCluster
from gmx.math.common import cluster, fitplane, convertSeconds
from gmx.math.hbond import hbondmolgrp
from gmx.math.rmsd import rmsdMol
from gmx.math.image import array2image, writeImage

from gmx.other.input_func import convert_input_type

from custom.general import listFiles
from custom.progressbar import ProgressBar

class TrajAnalysis(object):
	def __init__(self):
		self.trajprof = {}
		self.trajcoord = []
		self.trajn = 0
		self.atomn = 0
		self.moln = 0
		self.rmsdMatrix = np.array([])
		self.weightFactor = np.array([])

		self._savePath = "."
		self._loadPath = "."
		self._workPath = "."

		self.__rmsdMatrixCalc = False
		self.__outWeightFactor = False

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

	def getWorkPath(self):
		return self._workPath
	def setWorkPath(self, path):
		if os.path.isdir(path):
			self._workPath = path
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
		self.rmsdMatrix = np.array([])
		self.weightFactor = np.array([])

		# Reset Flags
		self.__rmsdMatrixCalc = False
		self.__outWeightFactor = False

	def useOutWeightFactor(self):
		return self.__outWeightFactor

	def checkRmsdMatrixCalc(self):
		return self.__rmsdMatrixCalc

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

	def loadWeightFactor(self):
		rightFactor = False
		file = listFiles(None, (".txt", ".csv", ".log"), self.getWorkPath())
		with open(file, "r") as f:
			factorLine = f.readline().strip().split()
			if len(factors) != self.atomn:
				rightFactor = False
			else:
				rightFactor = True
				for i in factorLine:
					try:
						float(i)
					except Exception as e:
						rightFactor = False
						break

		if rightFactor is True:
			factors = [float(i) for i in factorLine]
			self.weightFactor = np.array(factors)
			self.__outWeightFactor = True
		else:
			self.__outWeightFactor = False

	def setWeightFactor(self):
		rightFactor = False
		while (rightFactor is False):
			factorLine = input("Enter weighting factors for this trajectory: ")
			factors = factorLine.strip().split()
			if len(factors) != self.atomn:
				rightFactor = False
			else:
				rightFactor = True
				for i in factors:
					try:
						float(i)
					except Exception as e:
						rightFactor = False
						break
		factors = [float(i) for i in factors]
		self.weightFactor = np.array(factors)

		self.__outWeightFactor = True

	def loadRmsdMatrix(self):
		matrixFile = listFiles(None, ".csv", self.getWorkPath())

		matrix = []
		with open(matrixFile) as f:
			rightFormat = True
			columnN = -1
			for line in f.readlines():
				nums = []
				initline = line.strip().split(";")

				# Check columns number
				if columnN < 0:
					columnN = len(initline)
				elif columnN != len(initline):
					rightFormat = False
					break

				# Check whether each term in line is a real number
				for i in initline:
					try:
						n = float(i)
					except Exception as e:
						rightFormat = False
						break
					else:
						nums.append(n)
				if rightFormat is True:
					matrix.append(nums)
				else:
					break

		if rightFormat is True:
			matrix = np.array(matrix)
			if checkRmsdMatrix(matrix):
				self.rmsdMatrix = matrix
				self.trajn = matrix.shape[0]
				self.__rmsdMatrixCalc = True

				upperHalf = upperRmsdMatrix(self.rmsdMatrix)
				print("[Info] Minimum RMSD Value = %8.3f." % np.min(upperHalf))
				print("[Info] Maximum RMSD Value = %8.3f." % np.max(upperHalf))
				print("[Info] Average RMSD Value = %8.3f." % np.average(upperHalf))
			else:
				print("[Error] Wrong matrix file format. Maybe something is missing in the given matrix file. Please check your original matrix file.")
		else:
			print("[Error] Incompatible matrix file format.")

	def listCommand(self):
		cmdList = [x for x in dir(self) if not x.startswith("_") and callable(getattr(self, x))]
		cmdList.sort()
		for i in range(len(cmdList)):
			print(str(i + 1) + ": " + cmdList[i])

	def command(self):
		def exit():
			exitConfirm = input("Exit?(Y/N)").strip().upper()
			while exitConfirm not in ("Y", "N"):
				exitConfirm = input("Exit?(Y/N)").strip().upper()
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
			elif inputCmdList[0] == "command":
				print("[Info] Already in console mode.")
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
		ManageLoad = ProgressBar(self.trajn)
		ManageLoad.start("[Info] Start reading files.")
		for i in range(self.trajstartn, self.trajendn + 1):
			mols = pdb2molmatrix(framePath(i))
			self.trajcoord.append(np.vstack([mol.get("coordinate") for mol in mols.values()]))
			ManageLoad.forward(1)
		ManageLoad.end("[Info] All File have loaded successfully.")

		self.rmsdMatrix = np.zeros((self.trajn, self.trajn))
		self.weightFactor = np.zeros((self.trajn,))

	def writeLog(self, name, log, suffix = ".txt"):
		if not isinstance(name, str) or not (isinstance(log, (tuple, list)) and all([isinstance(l, str) for l in log])):
			raise TypeError("[Error] Wrong parameters. Function writeLog expects the first parameter to be a string and the second one to be a string list (or tuple).")
		else:
			# yyyymmdd-HHMMSS
			timeFormat = time.strftime("%Y%m%d_%H%M%S")
			fileName = name + "_" + timeFormat + suffix
			with open(os.path.join(self._savePath, fileName), "w") as f:
				f.writelines(line + '\n' for line in log)

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
			atomB = [atomB]
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
			atomid = [atomid]
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
			atomid = [atomid]
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
			molid = [molid]
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

	def dist2plane_mol_sign(self, molid, plane, reference_coord = (0., 0., 0.)):
		# Check
		if not self.checkMolID(molid) or not self.checkAtomIDList(plane):
			raise ValueError("Wrong parameters. Function dist2planeMol expects the first parameter to be the selected molecular ID, the second one to be an ID list (or tuple) for atom(s) in the reference plane and the last one to be the selected atom ID that should be in the positive site.")
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
			# Calculate reference distance value to determine the sign
			reference_dist = (plnormal.dot(reference_coord) + plps[3])/math.sqrt(plnormal.dot(plnormal))
			pldist = (plnormal.dot(center) + plps[3])/math.sqrt(plnormal.dot(plnormal))
			# Different sign
			if reference_dist * pldist < 0:
				pldist = - pldist
			log.append('{0:>6}'.format(framen) + " | " + '{0:>8.3f}'.format(pldist))

			self.writeLog("dist2plane_mol_sign", log)

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

	def concRmsdMatrix(self):
		# Check shape of RMSD matrix
		if self.rmsdMatrix.ndim == 2 and self.rmsdMatrix.shape[0] == self.trajn and self.rmsdMatrix.shape[1] == self.trajn:
			pass
		else:
			self.rmsdMatrix = np.zeros((self.trajn, self.trajn))

		# Construct RMSD Matrix

		# change weight factor
		if (self.useOutWeightFactor() is True):
			weightFactor = self.weightFactor
		else:
			weightFactor = np.ones(self.atomn)

		# Start Calculation
		startTime = time.time()
		ManageRMSD = ProgressBar(self.trajn - 1)
		ManageRMSD.start("[Info] Start calculating RMSD matrix.")
		for i in range(self.trajn - 1):
			for j in range(i + 1, self.trajn):
				self.rmsdMatrix[i][j] = rmsdMol(self.trajcoord[i], self.trajcoord[j], weightFactor)
				self.rmsdMatrix[j][i] = self.rmsdMatrix[i][j]
			ManageRMSD.forward(1)
		endTime = time.time()
		ManageRMSD.end("[Info] RMSD matrix calculation normally ends.")
		print("[Info] Total time usage = %s." % convertSeconds(endTime - startTime))

		upperHalf = upperRmsdMatrix(self.rmsdMatrix)
		print("[Info] Minimum RMSD Value = %8.3f." % np.min(upperHalf))
		print("[Info] Maximum RMSD Value = %8.3f." % np.max(upperHalf))
		print("[Info] Average RMSD Value = %8.3f." % np.average(upperHalf))

		rmsdLog = [";".join(['{0:>6.3f}'.format(i) for i in r]) for r in self.rmsdMatrix]

		self.writeLog("rmsd_matrix", rmsdLog, ".csv")

		self.__rmsdMatrixCalc = True

		# RMSD matrix image
		writeImage(
			array2image(self.rmsdMatrix, (self.trajn, self.trajn), colorMap = True),
			"rmsdMatrix",
			self.getSavePath()
		)

	def clusterSingleLinkage(self, cutoff):
		if self.checkRmsdMatrixCalc() is False:
			self.concRmsdMatrix()

		cindex, clust = renumCluster(
			linkStructure(
				self.trajn,
				cutoff,
				reshapeRmsdMatrix(self.rmsdMatrix)
			)
		)

		infoLine = "Method: Single Linkage, Cutoff: %.3f nm" % (cutoff)
		titleLine = '{0:>6}'.format("Groupn") + " | " + '{0:>24}'.format("Group members (Framen)")
		log = [infoLine, titleLine]
		groupn = 0
		for group in clust:
			groupn += 1
			log.append('{0:>6}'.format(groupn) + " | " + ', '.join([str(g) for g in group]))

		self.writeLog("cluster_single_linkage", log)

	def hbondgrp(self, moltypeA, moltypeB, lpdir = False):
		# Check
		profgrp = cluster(self.trajprof)
		if moltypeA not in profgrp.keys() or moltypeB not in profgrp.keys():
			raise ValueError("Wrong parameters. Function hbondgrp expects at least two parameters, both of them should be three-letter tags of corresponding molecules.")
		# Title
		frametitle = '{0:>14}'.format("donor") + " | " + '{0:>14}'.format("hydro") + " | " + '{0:>14}'.format("acceptor") + " | " + '{0:>9}'.format("distance") + " | " + '{0:>6}'.format("angle")
		log = []

		# Start Calculation
		startTime = time.time()
		ManageHbond = ProgressBar(self.trajn)
		ManageHbond.start("[Info] Start identifying hydrogen bond.")
		# Start Counting
		# Counting at the end of the loop
		framen = self.trajstartn
		frame_hbondns, frame_hbondmols = {}, {}
		for framecoord in self.trajcoord:
			hbondn, orglog, molids = hbondmolgrp(self.trajprof, framecoord, moltypeA, moltypeB, lpdir)
			frame_hbondn = sum(hbondn)
			frame_hbondns.update({framen: frame_hbondn})
			if frame_hbondn > 0:
				log.append("<frame %d>" % framen)
				log.append(frametitle)
				log.extend(orglog)
				frame_hbondmols.update({framen: molids})
			ManageHbond.forward(1)
			framen += 1
		# Calculation ends
		endTime = time.time()
		ManageHbond.end("[Info] Hygrogen bond identification normally ends.")
		print("[Info] Total time usage = %s." % convertSeconds(endTime - startTime))

		if log:
			self.writeLog("hbond_%s-%s" %(moltypeA.lower(), moltypeB.lower()), log)

		hbondnlog = ['{0:>6}'.format("frames") + " | " + '{0:>6}'.format("hbondn")]
		for i in range(self.trajstartn, self.trajendn + 1):
			hbondnlog.append('{0:>6}'.format(i) + " | " + '{0:>6}'.format(frame_hbondns.get(i)))
		self.writeLog("hbondsum_%s-%s" %(moltypeA.lower(), moltypeB.lower()), hbondnlog)

		timeFormat = time.strftime("%Y%m%d_%H%M%S")
		newFrameDir = os.path.join(self.getWorkPath(), "mols_%s-%s" %(moltypeA.lower(), moltypeB.lower())) + timeFormat
		if os.mkdir(newFrameDir):
			self.setWorkPath(newFrameDir)
			for j in frame_hbondmols.keys():
				self.writeNewFrame(j, frame_hbondmols.get(j))
			self.setWorkPath(os.path.abspath(os.path.join(newFrameDir), ".."))

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
		molmatrix2pdb(molmatrix, os.path.join(self.getWorkPath(), "newframe%d.pdb" % frameid))

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
			molmatrix2pdb(molmatrix, os.path.join(self.getWorkPath(), "%s.pdb" % filename))
		else:
			raise ValueError("Wrong parameters. Function writeGroups expects two parameters, the former one to be a tuple list in which each element should be a pack of frame ID and molecular ID list and the latter one to be the wanted file name without suffix.")
