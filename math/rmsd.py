import numpy as np

from gmx.structure.coordMatrix import checkCoordinateMatrix

def rmsd(refMatrix, modMatrix):
    # Check parameters
    status, errorInfo = checkCoordinateMatrix(refMatrix)
    if not status:
        raise ValueError(errorInfo)
    status, errorInfo = checkCoordinateMatrix(modMatrix)
    if not status:
        raise ValueError(errorInfo)


def fitting(refMatrix, modMatrix):
    # Check parameters
    status, errorInfo = checkCoordinateMatrix(refMatrix)
    if not status:
        raise ValueError(errorInfo)
    status, errorInfo = checkCoordinateMatrix(modMatrix)
    if not status:
        raise ValueError(errorInfo)
    refCenter = np.average(refMatrix, axis = 0)
    modCenter = np.average(modMatrix, axis = 0)
    # Translate the latter matrix to align centers of the two molecules
    translateMatrix = modMatrix + (refCenter - modCenter)
    # Rotate the latter matrix to fit with the least RMSD
    return newMatrix

def rmsdMol(molMatrixA, molMatrixB, fitting = True):
    pass
