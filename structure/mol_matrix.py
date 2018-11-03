import os.path

from gmx.io.pdb import pdb2mol_matrix, mol_matrix2pdb
from gmx.io.gro import gro2mol_matrix, mol_matrix2gro

def load_as_mol_matrix(frame_path):
	frame_name, frame_format = os.path.splitext(frame_path)
	mol_matrix = {}
	if frame_format == ".pdb":
		mol_matrix = pdb2mol_matrix(frame_name)
	elif frame_format == ".gro":
		mol_matrix = gro2mol_matrix(frame_name)
	return mol_matrix

def write_mol_matrix(mol_matrix, frame_path):
	frame_name, frame_format = os.path.splitext(frame_path)
	if frame_format == ".pdb":
		mol_matrix2pdb(mol_matrix, frame_name)
	elif frame_format == ".gro":
		mol_matrix2gro(mol_matrix, frame_name)
