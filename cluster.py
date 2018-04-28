# Import frequent used modules
import os
import gromacs as gmx
# Main class Cluster
class Cluster(object):
    def __init__(self):
        pass
    def load(self, dirpath, start, end, name, frametype = "pdb"):
		# Initialize
		self.trajProfile = {}
		self.trajCoordinate = []
		self.trajNum = end - start + 1
		# Load trajectory profile from the first frame
		initialFrame = pdb2molmatrix(os.path.join(dirpath, name + str(start) + "." + frametype))
		# Remove coordinate and velocity terms
		# Add atom id term for each atom, begin counting
		# Atom ids would be count from 1
		atomn = 1
		moln = 0
		for molID in initialFrame.keys():
			mol = initialFrame.get(molID)
			molAtomNum = len(mol.get("atom"))
			molProfile = {"name":mol.get("name"), "atom":mol.get("atom"), "atomid":[atomn+i for i in range(molAtomNum)]}
			atomn += molAtomNum
			moln += 1
			self.trajprof.update({molID:molProfile})
		# atomn equals to total counting of atoms plus initial 1
		self.atomNum = atomn - 1
		# moln equals to total molecule number
		self.molNum = moln
		# Only load coordinates data of each frame (including the first one) into a matrix (N*3)
		for i in range(start, end+1):
			molMatrix = pdb2molmatrix(os.path.join(dirpath, name + repr(i) + "." + frametype))
			self.trajcoord.append(np.vstack([mol.get("coordinate") for mol in molMatrix.values()]))
			print("[Info] File %s is loaded successfully." % (name + repr(i) + "." + frametype))
