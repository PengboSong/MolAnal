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
    refCenter = np.average(refMatrix*weightFactor, axis = 0)
    modCenter = np.average(modMatrix*weightFactor, axis = 0)
    # Translate two matrixes to align centers of the two molecules to the origin
    transRefMatrix = modMatrix - modCenter
    transModMatrix = modMatrix - modCenter
    # Rotate the latter matrix to fit with the least RMSD
    return newMatrix

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
        molNewMatrixB = fitting(molMatrixA, molMatrixB, weightFactor)
        rmsdValue = rmsd(molMatrixA, molNewMatrixB, weightFactor)
    else:
        rmsdValue = rmsd(molMatrixA, molMatrixB, weightFactor)
    return rmsdValue
