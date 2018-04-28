import numpy as np
from gmx.structure.coordMatrix import checkCoordinateMatrix
def fitting(refMatrix, modMatrix):
    # Check parameters
    status, errorInfo = checkCoordinateMatrix(refMatrix)
    if not status:
        raise ValueError(errorInfo)
    status, errorInfo = checkCoordinateMatrix(modMatrix)
    if not status:
        raise ValueError(errorInfo)



def rmsdMol(molMatrixA, molMatrixB):
    pass
