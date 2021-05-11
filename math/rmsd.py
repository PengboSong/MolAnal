# coding=utf-8

import math

import numpy as np
from gmx.other.data_type import GMXDataType


def weight(x, weight_factor=None):
    """Perform weighted average on x.
    
    Args:
        x: Data series to average. If x is a vector, x should have the same
           length with the weight factor. If x is a matrix, the first dimension
           of x should be the same as the weight factor length.
        weight_factor: Input weight factor.
    
    Returns:
        Weighted average result of x.
    """
    if isinstance(weight_factor, np.ndarray) and weight_factor.size == x.shape[0]:
        weight_factor /= np.sum(weight_factor)
        res = np.dot(weight_factor, x)
    else:
        res = np.mean(x, axis=0)
    return res


def rmsd(ref_mat, mol_mat, weight_factor=None):
    """Calculate RMSD value without fitting.

    Args:
        ref_mat: Referenced 2-dim coordinate matrix with 3 columns.
        mol_mat: Current 2-dim coordinate matrix with 3 columns.
        weight_factor: Weight factor contributed by each atoms. If weight
                       factor is zero, the atom will have no contributions to
                       the calculated RMSD score. If weight factor is not
                       specified, s simple average will be performed.
    
    Returns:
        Calculated RMSD value.
    """
    # Calculate length of random dimensions vector
    squared_length = np.sum(np.power(mol_mat - ref_mat, 2), axis=1)
    # RMSD = {1/N Σ(x-x0)^2 + (y-y0)^2 + (z-z0)^2}^(1/2)
    # N - Total number of atoms
    # x, y, z - X, Y, Z coordinates of atoms in current state
    # x0, y0, z0 - X, Y, Z coordinates of atoms in referenced state
    return math.sqrt(weight(squared_length, weight_factor))


def fitting(ref_mat, mol_mat, weight_factor):
    """Calculate RMSD value with fitting.

    Args:
        ref_mat: Referenced 2-dim coordinate matrix with 3 columns.
        mol_mat: Current 2-dim coordinate matrix with 3 columns.
        weight_factor: Weight factor contributed by each atoms. If weight
                       factor is zero, the atom will have no contributions to
                       the calculated RMSD score. If weight factor is not
                       specified, s simple average will be performed.
    
    Returns:
        Calculated RMSD value.
    """
    # Calculate vectors from origin to centers of molecules
    ref_center = weight(ref_mat, weight_factor)
    mol_center = weight(mol_mat, weight_factor)

    # Translate two matrixes to align centers of the two molecules to the origin
    trans_ref_matrix = ref_mat - ref_center
    trans_mod_matrix = mol_mat - mol_center

    # Rotate the latter matrix to fit with the least RMSD

    # Derive unitary rotation matrix
    NDIM = 3
    vectorX = trans_ref_matrix.transpose()
    vectorY = trans_mod_matrix.transpose()
    lenVector = vectorX.shape[1]

    rM = np.zeros((NDIM, NDIM), dtype=GMXDataType.REAL)

    for i in range(NDIM):
        for j in range(NDIM):
            for k in range(lenVector):
                rM[i][j] += weight_factor[k] * vectorY[i][k] * vectorX[j][k]

    sPDM = np.dot(rM.transpose(), rM)

    eigValues, eigVectors = np.linalg.eig(sPDM)

    bM = np.dot(rM, eigVectors) / np.sqrt(eigValues)
    aM = eigVectors.transpose()

    uM = np.zeros((NDIM, NDIM), dtype=GMXDataType.REAL)  # unitary rotation matrix
    # UX = Y
    # U'Y = X

    for i in range(NDIM):
        for j in range(NDIM):
            for k in range(NDIM):
                uM[i][j] += bM[k][j] * aM[k][i]

    new_matrix = np.dot(uM, vectorY).transpose() + ref_center
    fit_rmsd = rmsd(ref_mat, new_matrix, weight_factor)
    return new_matrix, fit_rmsd


def rmsd_mol(coord_mat1, coord_mat2, weight_factor, fitting=True):
    """Calculate RMSD value with given coordinate matrix and weight factor.

    Args:
        coord_mat1: Coordinate matrix for referenced state of molecule/structure.
        coord_mat2: Coordinate matrix for state of molecule/structure.
        fitting: If true, the two molecules/structures will be fitted before
                 calculating RMSD. Otherwise, calculate RMSD using original
                 coordinates.
    
    Returns:
        RMSD value between two molecules/structures.
    """
    # Check parameters
    weight_factor = np.asarray(weight_factor, dtype=GMXDataType.REAL).reshape(-1)
    if coord_mat1.atomn == coord_mat2.atomn and coord_mat1.atomn == weight_factor.size:
        if fitting:
            rot_coord_mat2, rmsd_value = fitting(coord_mat1, coord_mat1, weight_factor)
        else:
            rmsd_value = rmsd(coord_mat1, coord_mat2, weight_factor)
        return rmsd_value
    else:
        raise ValueError(
            "Lengths of reference matrix, modified matrix and weight factor vector are not equal.")
