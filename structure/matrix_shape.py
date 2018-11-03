import numpy as np

def is_coord_matrix(matrix):
	acceptableDataType = ["int32", "int64", "float64", "float128"]
	errorInfo = ""
	if isinstance(matrix, np.ndarray):
		if str(matrix.dtype) in acceptableDataType:
			if matrix.ndim == 2:
				if matrix.shape[1] == 3:
					return True, errorInfo
				else:
					errorInfo = "Unsupported data format. Only handle matrix with 3 columns."
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

def is_vector(matrix):
	acceptableDataType = ["int32", "int64", "float64", "float128"]
	errorInfo = ""
	if isinstance(matrix, np.ndarray):
		if str(matrix.dtype) in acceptableDataType:
			if matrix.ndim == 1:
				return True, errorInfo
			else:
				errorInfo = "Unsupported data format. Only handle 1-dim matrix(vector)."
				return False, errorInfo
		else:
			errorInfo = "Unsupported data format. Only handle matrix with number(integer or float) inside."
			return False, errorInfo
	else:
		errorInfo = "Unsupported data format. Only handle matrix of numpy.ndarray type."
		return False, errorInfo

def is_squared_matrix(matrix):
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
