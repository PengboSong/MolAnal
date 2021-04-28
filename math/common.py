# coding=utf-8

import math

import numpy as np


def boxpara(arrayxyz, extend=1.10):
    """Calculating solvation box parameters.

    Args:
            arrayxyz: A 2-dim array containing XYZ of all atoms in system.
            extend: Factor that controls how large the solvation box will be.
                    Should be no less than 1, otherwise the box is not large
                            enough to contain all atoms.

    Returns:
            Length of solvation box in X, Y and Z axes.
    """
    # Extend should be no less than 1
    extend = max(extend, 1.0)
    # Check shape of arrayxyz
    if arrayxyz.ndim != 2 or arrayxyz.shape[1] != 3:
        print('[WARNING] "boxpara" gets wrong paramters.'
              '"arrayxyz" should be a 2-dim matrix with 3 columns.')
        return 0., 0., 0.

    xmin, ymin, zmin = np.min(arrayxyz, axis=0)
    xmax, ymax, zmax = np.max(arrayxyz, axis=0)
    return (xmax-xmin) * (1+extend), (ymax-ymin) * (1+extend), (zmax-zmin) * (1+extend)


def fitplane(pts):
    """Fitting plane with given point coordinates.

    This function can be used to fitting a plane close to all points given.
    If only coordinates of 3 points given, fitplane will find a plane that
    all 3 points are just on the plane.
    If coordinates of more than 3 points given, firplane will find a plane
    that the angles between normal vector of the plane and vector from center
    to each points are close to 90 degrees (perpendicular).

    Args:
            pts: Coordinates of at least 3 points. Should be a 2-dim matrix with
                 3 columns being X, Y, Z coordinates.

    Returns:
            Fitting parameters of the plane: a, b, c, d.
            The Fitted plane can be expressed as equation "ax+by+cz+d=0"
    """
    # Use leastsq module from scipy to converge the least squared summation
    from scipy.optimize import leastsq

    def planeq(p, x, y, z):
        """Function definded for optimization of leastsq module.

        Args:
                p: Initial guess of plane normal vector. Should be a normal vector.
                x, y, z: X, Y, Z coordinates of fitting points.

        Returns:
                Distances from 
        """
        p = np.asarray(p, dtype="float64")
        pt = np.asarray([x, y, z], dtype="float64")
        # CAN NOT use np.dot(pt, pt) to replace (x**2 + y**2 + z**2) - They are NOT equal
        return (np.dot(p, pt))/np.sqrt((p.dot(p) * (x**2 + y**2 + z**2)))

    pts = np.asarray(pts, dtype="float64")
    if pts.ndim != 2 or pts.shape[0] < 3 or pts.shape[1] != 3:
        raise ValueError('"fitplane" gets wrong parameters.'
                         '"pts" should be a 2-dim matrix with 3 columns and at least 3 rows.')

    # Get the initial guess of normal vector from the first three points
    normal = np.cross(pts[1]-pts[0], pts[2]-pts[0])
    a, b, c = normal
    if len(pts) == 3:
        # When only 3 points given, plane equation can be determined precisely
        d = -(normal*pts[0]).sum()
        return a, b, c, d
    else:
        # When more than 3 points given, find a close plane
        # adjustPoints is a Nx3 matrix that each row represents a vector from center to a point
        adjustPoints = np.array(pts) - np.mean(np.array(pts), axis=0)
        # Transpose of adjustPoints is a 3xN matrix, three rows cooresponds to X, Y, Z, respectively
        x, y, z = adjustPoints.transpose()
        a, b, c = leastsq(planeq, normal, args=(x, y, z))[0]
        d = -np.array([a, b, c]).dot(np.mean(np.array(pts), axis=0))
        return a, b, c, d


def rotate(vector, axis, cosang, sinang):
    """Rotate a vector along rotation axis by angle specified by cos and sin.

    Args:
            vector: The vector to be rotated.
            axis: Vector of rotation axis. Should be a unit vector.
            cosang: Cosine of rotation angle.
            sinang: Sine of rotation angle.

    Returns:
            The vector after rotation, which has the same shape as input vector.
    """
    rotmatrix = rotateMatrix(axis, cosang, sinang)
    return np.dot(rotmatrix, vector)


def rotateMatrix(axis, cosang, sinang):
    """Construct a rotation matrix.

    Args:
            axis: Vector of rotation axis. Should be a unit vector.
            cosang: Cosine of rotation angle.
            sinang: Sine of rotation angle.

    Returns:
            A 3x3 rotation matrix, which should be a unitary matrix.
    """
    axis = axis / math.sqrt(axis.dot(axis))
    x, y, z = axis
    return (1 - cosang) * np.outer(axis, axis) + \
        sinang * np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]]) + \
        cosang * np.eye(3)


def convertSeconds(sec):
    """Convert seconds number to string listing hours, minutes and seconds.

    Args:
            sec: Time length in seconds. Should be a float.

    Returns:
            A string listing hours, minutes and seconds (correct to centiseconds).
    """
    timestr = ""
    h = int(sec // 3600)
    timestr += f"{h:d} hour" if h > 0 else ""
    lefts = sec % 3600
    m = int(lefts // 60)
    timestr += f"{m:d} min" if m > 0 else ""
    s = lefts % 60
    timestr += f"{s:.2f} sec"
    return timestr
