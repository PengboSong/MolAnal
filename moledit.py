import os, math, random
import numpy as np

from gmx.other.input_func import convert_input_type
from gmx.structure.mol_matrix import load_mol_matrix, write_mol_matrix

from gmx.math.common import fitplane, rotate

from custom.general import listFiles

class MolEdit(object):
	def __init__(self, frame):
		self.mols = load_mol_matrix(frame)
		self.parent_dir, self.frame_name = os.path.split(frame)
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
		if copy:
			frame_name = input("Please enter the name to save as:")
		else:
			frame_name = self.frame_name

		write_mol_matrix(self.mols, os.path.join(self.parent_dir, frame_name), True)
		print("[Info] File %s has been written successfully." % frame_name)
