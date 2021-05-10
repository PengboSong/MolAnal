# coding=utf-8

import math
import os
import random

import numpy as np
from custom.general import listFiles
from gmx.math.common import fitplane, rotate
from gmx.other.data_type import GMXDataType
from gmx.other.input_func import convert_input_type
from gmx.other.register_func import RegisterFunction
from gmx.structure.mol_matrix import MolMatrix


class MolEdit(RegisterFunction):
    PROMPT = "Molecules Editing for GROMACS"
    EXIT_PROMPT = "Press any key to exit..."
    PLACEHOLDER = ">>> "

    def __init__(self, frame=''):
        super().__init__()
        self.mols = MolMatrix()
        if frame:
            self.mols.from_file(frame)
        self.parent_dir = os.path.dirname(frame)
        self.frames = [os.path.basename(frame)]
    
    def __call__(self):
        """Initialize command-line interactive interface"""
        print(self.PROMPT)
        cmd = input(self.PLACEHOLDER).strip()
        lcmd = {"move":self.move, "moveto":self.moveto, "rotate":self.orrient, "random-rotate":self.random_orrient,"write":self.write}
        # Read commands from console
        while cmd != "exit":
            if cmd.split(" ")[0] not in lcmd.keys():
                print("[Info] Unknown command.")
            else:
                try:
                    lcmd.get(cmd.split(" ")[0])(*([convert_input_type(x) for x in cmd.split(" ")[1:]]))
                except Exception as e:
                    print("[Error] "+str(e)+".")
            cmd = input(self.PLACEHOLDER).strip()
        save_and_exit = input("Save modifications?(Y/N)").strip().upper()
        while save_and_exit not in ("Y", "N"):
            save_and_exit = input("Save modifications?(Y/N)").strip().upper()
        if save_and_exit == "Y":
            self.write()

    def load_mol(self, frame=''):
        """Load another molecule/frame file into current molecule matrix"""
        new_mol = MolMatrix()
        if frame:
            new_mol.from_file(frame)
        self.mols.merge(new_mol)
    
    def check_molid(self, molid):
        if molid > 0 and molid <= len(self.mols):
            return True
        else:
            return False
    
    def cut(self, centermol = 1, shape = "sphere", parameter = (1.0,)):
        """Cut molecules near central molecule within a custom shape."""
        # Check
        if self.check_molid(centermol):
            raise ValueError("Target molecule %d to cut not found." % centermol)
        # Calculate coordinate of target molecular center
        center = self.mols[centermol - 1].center()
        # Cut a sphere around target molecular center
        if shape == "sphere":
            inside = MolMatrix()
            outside = MolMatrix()
            radius = parameter[0]
            # Scan all molecules
            for mol in self.mols:
                if mol.i == centermol:
                    continue
                else:
                    molcenter = mol.center()
                    # Put a molecule "inside" when its center is within the sphere
                    if np.linalg.norm(molcenter - center) < radius:
                        inside.merge(mol)
                    else:
                        outside.merge(mol)
            return inside, outside
        else:
            raise ValueError("Unsupported shape name: %s." % shape)

    def move(self, movemol = 1, vector = [0., 0., 0.]):
        """Move molecule by a specified vector"""
        # Check
        vector = np.asarray(vector, dtype=GMXDataType.REAL).reshape(-1)
        assert vector.size == 3, "Moving vector should be like [vx, vy, vz]."
        if not self.check_molid(movemol):
            raise ValueError("Target molecule %d to move not found." % movemol)
        # Move
        self.mols[movemol].move(*vector)

    def moveto(self, movemol = 1, point = [0., 0., 0.]):
        """Move molecule center to a specified point"""
        # Check
        point = np.asarray(point, dtype=GMXDataType.REAL).reshape(-1)
        assert point.size == 3, "Moving target point should be like [x, y, z]."
        if not self.check_molid(movemol):
            raise ValueError("Target molecule %d to move not found." % movemol)
        # Move target molecule to destination
        # The molecule goes to where its center overlap the point
        movevec = point - self.mols[movemol].center()
        self.mols[movemol].move(*movevec)
    
    def orrient(self, rotmol = 1, axis = [0., 0., 1.]):
        # Check
        point = np.asarray(point, dtype=GMXDataType.REAL).reshape(-1)
        assert point.size == 3, "Moving target point should be like [x, y, z]."
        if not self.check_molid(rotmol):
            raise ValueError("Target molecule %d to move not found." % rotmol)
        # Normalize
        vector = vector/math.sqrt(vector.dot(vector))
        # Coordinate matrix of target molecule
        molcoord = self.mols.get(movemol).get("coordinate")
        center = np.average(molcoord, axis=0)
        # Get the normal vector determined by molecular coordinates
        plps = fitplane([v for v in molcoord])
        normal = np.array(plps[0:3])
        normal = normal/math.sqrt(normal.dot(normal))
        # Rotation axis
        axis = np.cross(normal, vector)
        axis = axis/math.sqrt(axis.dot(axis))
        # Rotation angle
        cosang = np.dot(normal, vector)
        sinang = math.sqrt(1-cosang**2)
        # Rotate
        newcoord = []
        for i in range(len(molcoord)):
            newcoord.append(rotate(molcoord[i]-center, axis, cosang, sinang) + center)
        self.mols.get(movemol).update({"coordinate":np.array(newcoord)})

    def random_orrient(self, rotmol):
        self.orrient(rotmol=rotmol, axis=np.random.random(3, dtype=GMXDataType.REAL))

    def write(self, copy = False):
        if copy:
            frame_name = input("Please enter the name to save as:")
        else:
            frame_name = self.frame_name

        write_mol_matrix(self.mols, os.path.join(self.parent_dir, frame_name), True)
        print("[Info] File %s has been written successfully." % frame_name)
