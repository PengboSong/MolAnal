import numpy as np

def checkRmsdMatrix(matrix):
    acceptableDataType = ["int32", "int64", "float64", "float128"]
    errorInfo = ""
    if isinstance(matrix, np.ndarray):
        if str(matrix.dtype) in acceptableDataType:
            if matrix.ndim == 2:
                if matrix.shape[0] == matrix.shape[1]:
                    return True, errorInfo
                else:
                    errorInfo = "Unsupported data format. Only handle symmetric matrix."
                    return False, errorInfo
            else:
                errorInfo = "Unsupported data format. Only handle 2-dim matrix."
                return False, errorInfo
        else:
            errorInfo = "Unsupported data format. Only handle matrix with number(integer or float) inside."
            return False, errorInfo
    else:
        errorInfo = "Unsupported data format. Only handle matrix of numpy.ndarray type."
        return False, errorInfo

def upperRmsdMatrix(matrix):
    sequence = []
    if checkRmsdMatrix(matrix):
        for i in range(matrix.shape[0]):
            sequence.extend(matrix[i][i + 1:].tolist())
    return np.array(sequence)

def reshapeRmsdMatrix(matrix):
    rmsdSequence = []
    if checkRmsdMatrix(matrix):
        N = matrix.shape[0]
        for i in range(N):
            for j in range(i + 1, N):
                rmsdSequence.append({"x":i, "y":j, "dist":matrix[i][j]})
        rmsdSequence.sort(key = lambda obj: obj.get("dist"))
        return rmsdSequence
    else:
        raise ValueError("Invalid format RMSD matrix received.")
        return None
