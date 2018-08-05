import math
import numpy as np

from gmx.structure.coordMatrix import checkCoordinateMatrix
from gmx.structure.coordMatrix import checkLineMatrix

def rmsd(refMatrix, modMatrix, weightFactor):
    # lambda function calculate length of random dimensions vector
    lengthSquared = np.array(list(map(lambda vector: sum([val**2 for val in vector]), modMatrix - refMatrix)))
    rmsdValue = math.sqrt(np.average(lengthSquared * weightFactor))
    return rmsdValue

def fitting(refMatrix, modMatrix, weightFactor):

    # Calculate vectors from origin to centers of molecules
    broadWeightFactor = np.array([weightFactor, weightFactor, weightFactor]).transpose()
    refCenter = np.average(refMatrix * broadWeightFactor, axis = 0)
    modCenter = np.average(modMatrix * broadWeightFactor, axis = 0)

    # Translate two matrixes to align centers of the two molecules to the origin
    transRefMatrix = refMatrix - refCenter
    transModMatrix = modMatrix - modCenter

    # Rotate the latter matrix to fit with the least RMSD

    # Derive unitary rotation matrix
    ndim = 3
    vectorX = transRefMatrix.transpose()
    vectorY = transModMatrix.transpose()
    lenVector = vectorX.shape[1]

    rM = np.zeros((ndim, ndim))

    for i in range(ndim):
    	for j in range(ndim):
    		for k in range(lenVector):
    			rM[i][j] += weightFactor[k] * vectorY[i][k] * vectorX[j][k]

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

    newMatrix = np.dot(uM.transpose(), vectorY).transpose() + refCenter
    fittingRMSD = rmsd(refMatrix, newMatrix, weightFactor)
    return newMatrix, fittingRMSD

def rmsdMol(molMatrixA, molMatrixB, weightFactor, fittingFlag = True):
    # Check parameters
    status, errorInfo = checkCoordinateMatrix(molMatrixA)
    if not status:
        raise ValueError(errorInfo)
    status, errorInfo = checkCoordinateMatrix(molMatrixB)
    if not status:
        raise ValueError(errorInfo)
    status, errorInfo = checkLineMatrix(weightFactor)
    if not status:
        raise ValueError(errorInfo)
    if not (molMatrixA.shape[0] == molMatrixB.shape[0] == len(weightFactor)):
        raise ValueError("Lengths of reference matrix, modified matrix and weight factor vector are not equal.")
    if fittingFlag is True:
        molNewMatrixB, rmsdValue = fitting(molMatrixA, molMatrixB, weightFactor)
    else:
        rmsdValue = rmsd(molMatrixA, molMatrixB, weightFactor)
    return rmsdValue
