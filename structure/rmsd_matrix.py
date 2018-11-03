import numpy as np

from gmx.structure.matrix_shape import is_squared_matrix

def upper_rmsd_matrix(matrix):
    sequence = []
    if is_squared_matrix(matrix):
        for i in range(matrix.shape[0]):
            sequence.extend(matrix[i][i + 1:].tolist())
    return np.array(sequence)

def reshape_rmsd_matrix(matrix):
    rmsdSequence = []
    if is_squared_matrix(matrix):
        N = matrix.shape[0]
        for i in range(N):
            for j in range(i + 1, N):
                rmsdSequence.append({"x":i, "y":j, "dist":matrix[i][j]})
        rmsdSequence.sort(key = lambda obj: obj.get("dist"))
        return rmsdSequence
    else:
        raise ValueError("Invalid format RMSD matrix received.")
        return None
