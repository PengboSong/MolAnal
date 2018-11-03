import math
import numpy as np

from gmx.structure.coord_matrix import is_coord_matrix, is_vector

def rmsd(ref_mat, mol_mat, weight_factor):
	# lambda function calculate length of random dimensions vector
	length_squared = np.array(list(map(lambda vector: sum([val**2 for val in vector]), mol_mat - ref_mat)))
	rmsd_value = math.sqrt(np.average(length_squared * weight_factor))
	return rmsd_value

def fitting(ref_mat, mol_mat, weight_factor):

	# Calculate vectors from origin to centers of molecules
	broad_weight_factor = np.array([weight_factor, weight_factor, weight_factor]).transpose()
	ref_center = np.average(ref_mat * broad_weight_factor, axis = 0)
	mod_center = np.average(mol_mat * broad_weight_factor, axis = 0)

	# Translate two matrixes to align centers of the two molecules to the origin
	trans_ref_matrix = ref_mat - ref_center
	trans_mod_matrix = mol_mat - mod_center

	# Rotate the latter matrix to fit with the least RMSD

	# Derive unitary rotation matrix
	ndim = 3
	vectorX = trans_ref_matrix.transpose()
	vectorY = trans_mod_matrix.transpose()
	lenVector = vectorX.shape[1]

	rM = np.zeros((ndim, ndim))

	for i in range(ndim):
		for j in range(ndim):
			for k in range(lenVector):
				rM[i][j] += weight_factor[k] * vectorY[i][k] * vectorX[j][k]

	sPDM = np.dot(rM.transpose(), rM)

	eigValues, eigVectors = np.linalg.eig(sPDM)

	bVectors = []
	for l in range(ndim):
		bVectors.append(np.dot(rM, eigVectors[:, l]) / math.sqrt(eigValues[l]))
	aM = eigVectors.transpose()
	bM = np.array(bVectors)

	uM = np.zeros((ndim, ndim)) # unitary rotation matrix
	# UX = Y
	# U'Y = X

	for i in range(ndim):
		for j in range(ndim):
			for k in range(ndim):
				uM[i][j] += bM[k][i] * aM[k][j]

	new_matrix = np.dot(uM.transpose(), vectorY).transpose() + ref_center
	fit_rmsd = rmsd(ref_mat, new_matrix, weight_factor)
	return new_matrix, fit_rmsd

def rmsd_mol(molMatrixA, molMatrixB, weight_factor, fittingFlag = True):
	# Check parameters
	status, errorInfo = is_coord_matrix(molMatrixA)
	if not status:
		raise ValueError(errorInfo)
	status, errorInfo = is_coord_matrix(molMatrixB)
	if not status:
		raise ValueError(errorInfo)
	status, errorInfo = is_vector(weight_factor)
	if not status:
		raise ValueError(errorInfo)
	if not (molMatrixA.shape[0] == molMatrixB.shape[0] == len(weight_factor)):
		raise ValueError("Lengths of reference matrix, modified matrix and weight factor vector are not equal.")
	if fittingFlag is True:
		molNewMatrixB, rmsd_value = fitting(molMatrixA, molMatrixB, weight_factor)
	else:
		rmsd_value = rmsd(molMatrixA, molMatrixB, weight_factor)
	return rmsd_value
