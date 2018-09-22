import math
import numpy as np

# Attention: extend should be greater than 1
def boxpara(arrayxyz, extend = 1.10):
	if isinstance(arrayxyz, np.ndarray) and arrayxyz.ndim == 2 and arrayxyz.shape[1] == 3 and isinstance(extend, (int, float)) and extend > 1:
		xmin, ymin, zmin = np.min(arrayxyz, axis = 0)
		xmax, ymax, zmax = np.max(arrayxyz, axis = 0)
		return (xmax-xmin) * (1+extend), (ymax-ymin) * (1+extend), (zmax-zmin) * (1+extend)
	else:
		raise ValueError("Get wrong paramters. \"arrayxyz\" should be a two-dim numpy.ndarray with three columns and \"extend\" should be greater than 1.")

# Group molecules
def cluster(mols):
	groups = {}
	for molid in range(1, len(mols)+1):
		moltype = mols.get(molid).get("name")
		# New molecular group
		if not groups:
			groups.update({moltype:[molid]})
		# Attribute the molecule to an existing group
		else:
			if moltype in groups.keys():
				groups.get(moltype).append(molid)
			else:
				groups.update({moltype:[molid]})
	return groups

# Function fitplane should be used to give the equation of a plane that determined precisely by three points or "close" to more than three points
# "close" here means that each cosine of angle between normal vector and the vector from the center of cluster of points to one point should be as close to 90 degrees (perpendicular) as possible
# Function fitplane return four parameters: a, b, c, d. The equation of wanted plane is "ax+by+cz+d=0"
def fitplane(pts):
	# Use leastsq module from scipy to converge the least squared summation
    from scipy.optimize import leastsq
	# Function planeq is used for leastsq
	# p should be a normal vector(to initialize) or objects can be converted to a vector
    def planeq(p, x, y, z):
        if isinstance(p, np.ndarray) and p.size == 3 and p.ndim == 1:
            pass
        elif isinstance(p, (tuple, list)) and len(p) == 3 and all([isinstance(v, (int, float)) for v in p]):
            p = np.array(p)
        else:
            raise TypeError("Unsupported parameter given for plane equation.")
        pt = np.array([x, y, z])
		# CAN NOT use np.dot(pt, pt) to replace (x**2 + y**2 + z**2) - They are NOT equal
        return (np.dot(p, pt))/np.sqrt((p.dot(p) * (x**2 + y**2 + z**2)))
    if isinstance(pts, (tuple, list)) and len(pts) > 2 and all([isinstance(pt, np.ndarray) and pt.size == 3 and pt.ndim == 1 for pt in pts]):
		# Get the normal vector from the first three points
        normal = np.cross(pts[1]-pts[0], pts[2]-pts[0])
        a, b, c = normal
        if len(pts) == 3:
			# When only three points are given, return the precise equation of the plane
            d = -(normal*pts[0]).sum()
            return a, b, c, d
        else:
			# When more than three points are given, just return equation of a "close" plane
			# adjustPoints is a N*3 matrix that each row equals to a vector from the center of N points to the point
            adjustPoints = np.array(pts) - np.average(np.array(pts), axis=0)
			# Transpose of adjustPoints is a 3*N matrix, three rows cooresponds to x, y, z, respectively
            x, y, z = adjustPoints.transpose()
            a, b, c = leastsq(planeq, normal, args=(x, y, z))[0]
            d = -np.array([a, b, c]).dot(np.average(np.array(pts), axis=0))
            return a, b, c, d
    else:
        raise ValueError("Unsupported parameter given for plane function in function fitplane from module gromacs.")

# Attention: axis vector must be a unit vector
# OR, the rotation matrix cannot be a unitary matrix
def rotate(vector, axis, cosang, sinang):
	rotmatrix = rotateMatrix(axis, cosang, sinang)
	return np.dot(rotmatrix, vector)

def rotateMatrix(axis, cosang, sinang):
	axis = axis / math.sqrt(axis.dot(axis))
	x, y, z = axis
	return (1-cosang)*np.outer(axis, axis) + sinang*np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]]) + cosang*np.eye(3)

def convertSeconds(timeLength):
	timeLength = round(timeLength)
	if (timeLength < 60):
		return str(timeLength) + " sec"
	elif (timeLength < 3600):
		minValue = timeLength // 60
		secValue = timeLength % 60
		return str(minValue) + " min " + str(secValue) + " sec"
	else:
		hourValue = timeLength // 3600
		minValue = (timeLength % 3600) // 60
		secValue = (timeLength % 3600) % 60
		return str(hourValue) + " hour " + str(minValue) + " min " + str(secValue) + " sec"
