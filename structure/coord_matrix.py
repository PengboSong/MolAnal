import os.path

from gmx.io.pdb import pdb2coord_matrix
from gmx.io.gro import gro2coord_matrix

def load_coord_matrix(frame_path):
	frame_name, frame_format = os.path.splitext(frame_path)
	if frame_format == ".pdb":
		return pdb2coord_matrix(frame_path)
	elif frame_format == ".gro":
		return gro2coord_matrix(frame_path)
	else:
		return None
