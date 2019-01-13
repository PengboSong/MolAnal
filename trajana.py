import os, time, math
import numpy as np

from gmx.structure.coord_matrix import load_coord_matrix
from gmx.structure.mol_matrix import load_mol_matrix, write_mol_matrix
from gmx.structure.matrix_shape import is_squared_matrix
from gmx.structure.rmsd_matrix import upper_rmsd_matrix, reshape_rmsd_matrix

from gmx.math.clust import link_structure, renum_cluster
from gmx.math.common import cluster, fitplane, convertSeconds
from gmx.math.hbond import hbondmolgrp
from gmx.math.rmsd import rmsd_mol
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
		self.rmsd_mat = np.array([])
		self.weight = np.array([])

		self._save_path = "."
		self._load_path = "."
		self._work_path = "."

		self.__calc_rmsd_mat_flag = False
		self.__set_weight_flag = False

	def get_save_path(self):
		return self._save_path
	def set_save_path(self, path):
		if os.path.isdir(path):
			self._save_path = path
		else:
			print("[Error] Invalid path. The given path should be a existing directory.")

	def get_load_path(self):
		return self._load_path
	def set_load_path(self, path):
		if os.path.isdir(path):
			self._load_path = path
		else:
			print("[Error] Invalid path. The given path should be a existing directory.")

	def get_work_path(self):
		return self._work_path
	def set_work_path(self, path):
		if os.path.isdir(path):
			self._work_path = path
		else:
			print("[Error] Invalid path. The given path should be a existing directory.")

	def clear_all(self):
		self.trajprof = {}
		self.trajcoord = []
		self.trajn = 0
		self.trajstartn = 0
		self.trajendn = 0
		self.atomn = 0
		self.moln = 0
		self.rmsd_mat = np.array([])
		self.weight = np.array([])

		# Reset Flags
		self.__calc_rmsd_mat_flag = False
		self.__set_weight_flag = False

	def set_weight_status(self):
		return self.__set_weight_flag

	def calc_rmsd_mat_status(self):
		return self.__calc_rmsd_mat_flag

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

	def load_weight_factor(self):
		right_weight_flag = False
		file = listFiles(None, (".txt", ".csv", ".log"), self.get_work_path())
		with open(file, "r") as f:
			factorLine = f.readline().strip().split()
			if len(factors) != self.atomn:
				right_weight_flag = False
			else:
				right_weight_flag = True
				for i in factorLine:
					try:
						float(i)
					except Exception as e:
						right_weight_flag = False
						break

		if right_weight_flag is True:
			factors = [float(i) for i in factorLine]
			self.weight = np.array(factors)
			self.__set_weight_flag = True
		else:
			self.__set_weight_flag = False

	def set_weight_factor(self):
		right_weight_flag = False
		while (right_weight_flag is False):
			factorLine = input("Enter weighting factors for this trajectory: ")
			factors = factorLine.strip().split()
			if len(factors) != self.atomn:
				right_weight_flag = False
			else:
				right_weight_flag = True
				for i in factors:
					try:
						float(i)
					except Exception as e:
						right_weight_flag = False
						break
		factors = [float(i) for i in factors]
		self.weight = np.array(factors)

		self.__set_weight_flag = True

	def load_rmsd_matrix(self):
		matrixFile = listFiles(None, ".csv", self.get_work_path())

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
			if is_squared_matrix(matrix):
				self.rmsd_mat = matrix
				self.trajn = matrix.shape[0]
				self.__calc_rmsd_mat_flag = True

				upperHalf = upperRmsdMatrix(self.rmsd_mat)
				print("[Info] Minimum RMSD Value = %8.3f." % np.min(upperHalf))
				print("[Info] Maximum RMSD Value = %8.3f." % np.max(upperHalf))
				print("[Info] Average RMSD Value = %8.3f." % np.average(upperHalf))
			else:
				print("[Error] Wrong matrix file format. Maybe something is missing in the given matrix file. Please check your original matrix file.")
		else:
			print("[Error] Incompatible matrix file format.")

	def list_command(self):
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

	def load(self, start, end, prefix = "frame", suffix = "pdb"):
		# Initialize
		self.trajn = end - start + 1
		self.trajstartn = start
		self.trajendn = end
		# Load trajectory profile from the first frame
		first_frame = load_mol_matrix(os.path.join(self.get_load_path(), "%s%d.%s" % (prefix, start, suffix)))
		# Remove coordinate and velocity terms
		# Add atom id term for each atom, begin counting
		# Atom ids would be count from 1
		atomn = 1
		moln = 0
		molids = list(first_frame.keys())
		molids.sort()
		for molid in molids:
			mol = first_frame.get(molid)
			mol_atomn = len(mol.get("atom"))
			mol_prof = {
			"name":mol.get("name"),
			"atom":mol.get("atom"),
			"atomid":[atomn + i for i in range(mol_atomn)]
			}
			atomn += mol_atomn
			moln += 1
			self.trajprof.update({molid:mol_prof})
		# atomn equals to total counting of atoms minus initial 1
		self.atomn = atomn - 1
		# moln equals to total molecule number
		self.moln = moln
		# Only load coordinates data of each frame (including the first one) into a matrix (N*3)
		ManageLoad = ProgressBar(self.trajn)
		ManageLoad.start("[Info] Start reading files.")
		for i in range(self.trajstartn, self.trajendn + 1):
			self.trajcoord.append(
				load_coord_matrix(os.path.join(self.get_load_path(), "%s%d.%s" % (prefix, i, suffix)))
			)
			ManageLoad.forward(1)
		ManageLoad.end("[Info] All File have loaded successfully.")

		self.rmsd_mat = np.zeros((self.trajn, self.trajn))
		self.weight = np.zeros((self.trajn,))

	def mol_center(self, frameID, molID):
		atomIDList = self.trajprof.get(molID).get("atomid")
		framecoord = self.trajprof[frameID - self.trajstartn]
		center = np.average(framecoord[min(atomIDList) - 1 : max(atomIDList)], axis=0)
		return center
	def mol_center_with_coord(self, framecoord, molID):
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

		with open(os.path.join(self.get_save_path(), "distance-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

	def distance_mol(self, molid, atomid):
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
			center = self.mol_center_with_coord(framecoord, molid)
			dists = []
			for a in atomid:
				pair0 = center - framecoord[a - 1]
				findist = math.sqrt(pair0.dot(pair0))
				dists.append(findist)
			log.append('{0:>6}'.format(framen) + " | " + " | ".join(['{0:>11.3f}'.format(d) for d in dists]))

		with open(os.path.join(self.get_save_path(), "distance_mol-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

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

		with open(os.path.join(self.get_save_path(), "location-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

	def location_mol(self, molid):
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
				centers.extend(self.mol_center_with_coord(framecoord, m))
			centers = np.array(centers)
			sum += centers
			log.append('{0:>8}'.format(framen) + " | " + " | ".join(['{0:>6.3f}'.format(c) for c in centers]))
		log.append('{0:>8}'.format("Average") + " | " + " | ".join(['{0:>6.3f}'.format(c / (framen + 1)) for c in sum]))

		with open(os.path.join(self.get_save_path(), "location_mol-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

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

		with open(os.path.join(self.get_save_path(), "dist2plane-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

	def dist2plane_mol(self, molid, plane):
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
			center = self.mol_center_with_coord(framecoord, molid)
			# Function fitplane only accepts a package of at least three points' coordinate (better near a plane) in the form of either tuple or list, each element of the package should be a 1-dim 3-size numpy.ndarray object
			plps = fitplane([framecoord[x - 1] for x in plane])	# Short of plane parameters
			# Attention: plps should be a tuple with 4 elements, denoted as a, b, c, d. The plane equation should be ax+by+cz=0
			plnormal = np.array(plps[0:3])
			pldist = abs(plnormal.dot(center) + plps[3])/math.sqrt(plnormal.dot(plnormal))
			log.append('{0:>6}'.format(framen) + " | " + '{0:>8.3f}'.format(pldist))

		with open(os.path.join(self.get_save_path(), "dist2plane_mol-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

	def dist2plane_mol_signed(self, molid, plane, reference_coord = (0., 0., 0.)):
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
			center = self.mol_center_with_coord(framecoord, molid)
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

		with open(os.path.join(self.get_save_path(), "dist2plane_mol_sign-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

	def rmsd(self, framen):
		# Check
		if not self.checkFrameID(framen):
			raise ValueError("Wrong parameters. Function rmsd expects one parameter to be the selected frame ID.")

		# Title
		titleLine = '{0:>6}'.format("frames") + " | " + '{0:>8}'.format("RMSD")
		log = [titleLine]
		# Reference frame
		refMatrix = self.trajcoord[framen - self.trajstartn]
		weight_factor = np.ones(self.atomn)
		# Start Counting
		framen = self.trajstartn
		for i in range(self.trajn):
			framen += 1
			thisMatrix = self.trajcoord[i]
			rmsdv = rmsd_mol(refMatrix, thisMatrix, weight_factor)
			log.append('{0:>6}'.format(framen) + " | " + '{0:>8.3f}'.format(rmsdv))

		with open(os.path.join(self.get_save_path(), "rmsd_mol-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

	def conc_rmsd_matrix(self):
		# Check shape of RMSD matrix
		if self.rmsd_mat.ndim == 2 and self.rmsd_mat.shape[0] == self.trajn and self.rmsd_mat.shape[1] == self.trajn:
			pass
		else:
			self.rmsd_mat = np.zeros((self.trajn, self.trajn))

		# Construct RMSD Matrix

		# change weight factor
		if (self.set_weight_status() is True):
			weight_factor = self.weight
		else:
			weight_factor = np.ones(self.atomn)

		# Start Calculation
		startTime = time.time()
		ManageRMSD = ProgressBar(self.trajn - 1)
		ManageRMSD.start("[Info] Start calculating RMSD matrix.")
		for i in range(self.trajn - 1):
			for j in range(i + 1, self.trajn):
				self.rmsd_mat[i][j] = rmsd_mol(self.trajcoord[i], self.trajcoord[j], weight_factor)
				self.rmsd_mat[j][i] = self.rmsd_mat[i][j]
			ManageRMSD.forward(1)
		endTime = time.time()
		ManageRMSD.end("[Info] RMSD matrix calculation normally ends.")
		print("[Info] Total time usage = %s." % convertSeconds(endTime - startTime))

		upperHalf = upper_rmsd_matrix(self.rmsd_mat)
		print("[Info] Minimum RMSD Value = %8.3f." % np.min(upperHalf))
		print("[Info] Maximum RMSD Value = %8.3f." % np.max(upperHalf))
		print("[Info] Average RMSD Value = %8.3f." % np.average(upperHalf))

		with open(os.path.join(self.get_save_path(), "rmsd_matrix-%s.csv" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(";".join(['{0:>6.3f}'.format(i) for i in r]) + '\n' for r in self.rmsd_mat)

		self.__calc_rmsd_mat_flag = True

		# RMSD matrix image
		writeImage(
			array2image(self.rmsd_mat, (self.trajn, self.trajn), colorMap = True),
			"rmsd_matrix",
			self.get_save_path()
		)

	def cluster_single_linkage(self, cutoff):
		if self.calc_rmsd_mat_status() is False:
			self.concRmsdMatrix()

		cindex, clust = renum_cluster(
			link_structure(
				self.trajn,
				cutoff,
				reshape_rmsd_matrix(self.rmsd_mat)
			)
		)

		infoLine = "Method: Single Linkage, Cutoff: %.3f nm" % (cutoff)
		titleLine = '{0:>6}'.format("Groupn") + " | " + '{0:>24}'.format("Group members (Framen)")
		log = [infoLine, titleLine]
		groupn = 0
		for group in clust:
			groupn += 1
			log.append('{0:>6}'.format(groupn) + " | " + ', '.join([str(g) for g in group]))

		with open(os.path.join(self.get_save_path(), "cluster_single_linkage-%s.txt" % time.strftime("%Y%m%d_%H%M")), "w") as f:
			f.writelines(line + '\n' for line in log)

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
			with open(os.path.join(self.get_save_path(),
			"hbond_%s+%s-%s.txt" % (moltypeA.lower(), moltypeB.lower(), time.strftime("%Y%m%d_%H%M"))
			), "w") as f:
				f.writelines(line + '\n' for line in log)

		hbondnlog = ['{0:>6}'.format("frames") + " | " + '{0:>6}'.format("hbondn")]
		for i in range(self.trajstartn, self.trajendn + 1):
			hbondnlog.append('{0:>6}'.format(i) + " | " + '{0:>6}'.format(frame_hbondns.get(i)))
		with open(os.path.join(self.get_save_path(),
			"hbondsum_%s+%s-%s.txt" % (moltypeA.lower(), moltypeB.lower(), time.strftime("%Y%m%d_%H%M"))
			), "w") as f:
			f.writelines(line + '\n' for line in hbondnlog)

		new_frame_dir = os.path.join(self.get_save_path(), "mols_%s-%s" %(moltypeA.lower(), moltypeB.lower())) + time.strftime("%Y%m%d")
		if not os.path.isdir(new_frame_dir):
			os.mkdir(newFrameDir)

		org_save_path = self.get_save_path()
		self.set_save_path(new_frame_dir)
		for j in frame_hbondmols.keys():
			self.write_new_frame(j, frame_hbondmols.get(j))
		self.set_save_path(org_save_path)

	# Extract some molecules out of one frame
	def write_new_frame(self, frameid, molid, conect = False):
		molid = list(set(molid))
		molid.sort()
		if not self.checkFrameID(frameid) or not self.checkMolIDList(molid):
			raise ValueError("Wrong parameters. Function write_new_frame expects two parameters, the former one to be the selected frame ID and the latter one to be an ID list (or tuple) for the wanted molecule(s).")
		molmatrix = {}
		framecoord = self.trajcoord[frameid]
		moln = 0
		# Reconstruct the molecular matrix
		for i in molid:
			moln += 1
			mol_prof = self.trajprof.get(i)
			mol = {"name":mol_prof.get("name"), "atom":mol_prof.get("atom")}
			mol.update({"coordinate":np.vstack([framecoord[l-1] for l in mol_prof.get("atomid")])})
			molmatrix.update({moln:mol})
		write_mol_matrix(molmatrix, os.path.join(self.get_save_path(), "newframe%d.pdb" % frameid), conect)

	def write_groups(self, pairs, filename, conect = False):
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
					mol_prof = self.trajprof.get(i)
					mol = {"name":mol_prof.get("name"), "atom":mol_prof.get("atom")}
					mol.update({"coordinate":np.vstack([self.trajcoord[frameid][l-1] for l in mol_prof.get("atomid")])})
					molmatrix.update({moln:mol})
			write_mol_matrix(molmatrix, os.path.join(self.get_save_path(), "%s.pdb" % filename), conect)
		else:
			raise ValueError("Wrong parameters. Function write_groups expects two parameters, the former one to be a tuple list in which each element should be a pack of frame ID and molecular ID list and the latter one to be the wanted file name without suffix.")
