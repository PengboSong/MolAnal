# coding=utf-8

import math
import os

import numpy as np
from custom.general import listFiles

from gmx.math.common import fitplane
from gmx.other.command_line import ConsoleBase
from gmx.other.data_type import GMXDataType
from gmx.other.register_func import RegisterFunction
from gmx.structure.mol_matrix import MolMatrix


class MolEdit(ConsoleBase, RegisterFunction):
    PROMPT = "Molecules Editing for GROMACS"
    EXIT_PROMPT = "Press any key to exit..."
    PLACEHOLDER = ">>> "

    def __init__(self, frame=''):
        super().__init__()
        self.mols = MolMatrix()
        if frame:
            self.mols.from_file(frame)
        self.frames = [frame]

        self.register_("load_mol", self.load_mol, str)
        self.register_("move", self.move, int, list)
        self.register_("moveto", self.moveto, int, list)
        self.register_("orrient", self.orrient, int, list)
        self.register_("random_orrient", self.random_orrient, int)
        self.register_("save", self.save)
        self.register_("saveas", self.saveas, str)
        self.register_("help", self.help_)
        self.register_("exit", self.exit_)
    
    def __call__(self):
        """Initialize command-line interactive interface"""
        print(self.PROMPT)
        # Read commands from console
        EXIT_PATTERN = re.compile(r"[eE][xX][iI][tT]")
        while not self.exit_signal_:
            try:
                self.launch_(input(self.PLACEHOLDER))
            except Exception as err:
                print("[Error] {}".format(err))
        if self.askyn(prompt="Save modifications?"):
            self.save()

    def load_mol(self, frame=''):
        """Load another molecule/frame file into current molecule matrix.
        
        Args:
            frame(str): Path to molecule/frame file.
        """
        if os.path.isfile(frame):
            new_mol = MolMatrix()
            new_mol.from_file(frame)
            self.mols.merge(new_mol)
            self.frames.append(frame)
    
    def check_molid(self, molid):
        if molid > 0 and molid <= len(self.mols):
            return True
        else:
            return False
    
    def cut(self, centermol, shape="sphere", parameter=(1.0,)):
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

    def move(self, movemol, vector):
        """Move molecule by a specified vector.
        
        Args:
            movemol(int): Index of molecule to move (start from 1).
            vector(list[int*3]): Vector [vx, vy, vz] to move.

        Raise:
            ValueError: Molecule not found.
        """
        # Check
        vector = np.asarray(vector, dtype=GMXDataType.REAL).reshape(-1)
        assert vector.size == 3, "Moving vector should be like [vx, vy, vz]."
        if not self.check_molid(movemol):
            raise ValueError("Target molecule %d to move not found." % movemol)
        # Move
        self.mols[movemol].move(*vector)

    def moveto(self, movemol, point):
        """Move molecule center to a specified point.
        
        Args:
            movemol(int): Index of molecule to move (start from 1).
            point(list[int*3]): Destination coordinate [x, y, z] of molecule center.

        Raise:
            ValueError: Molecule not found.
        """
        # Check
        point = np.asarray(point, dtype=GMXDataType.REAL).reshape(-1)
        assert point.size == 3, "Moving target point should be like [x, y, z]."
        if not self.check_molid(movemol):
            raise ValueError("Target molecule %d to move not found." % movemol)
        # Move target molecule to destination
        # The molecule goes to where its center overlap the point
        movevec = point - self.mols[movemol].center()
        self.mols[movemol].move(*movevec)
    
    def orrient(self, rotmol, vector):
        """Rotate molecule to a specified direction.
        
        Args:
            rotmol(int): Index of molecule to rotate (start from 1).
            axis(list[int*3]): Vector [vx, vy, vz] as the rotation direction.

        Raise:
            ValueError: Molecule not found.
        """
        # Check
        axis = np.asarray(axis, dtype=GMXDataType.REAL).reshape(-1)
        assert axis.size == 3, "Rotation axis should be like [ax, ay, az]."
        if not self.check_molid(rotmol):
            raise ValueError("Target molecule %d to rotate not found." % rotmol)
        # Normalize
        axis /= np.linalg.norm(axis)
        # Get the normal vector determined by molecular coordinates
        plparms = fitplane(self.mols[rotmol].xyzs)[:3]
        normal = np.asarray(plparms, dtype=GMXDataType.REAL)
        normal /= np.linalg.norm(normal)
        # Rotation axis
        axis = np.cross(normal, vector)
        axis /= np.linglg.norm(axis)
        # Rotation angle
        cosang = np.dot(normal, vector)
        sinang = math.sqrt(1 - cosang**2)
        # Rotate
        self.mols[rotmol].rotate(axis, cosang, sinang)

    def random_orrient(self, rotmol):
        """Rotate molecule to a random direction.
        
        Args:
            rotmol(int): Index of molecule to rotate (start from 1).

        Raise:
            ValueError: Molecule not found.
        """
        self.orrient(rotmol=rotmol, axis=np.random.random(3, dtype=GMXDataType.REAL))
    
    def save(self):
        """Save to the first loaded molecule/frame file"""
        if self.askyn(prompt="Confirm overwrite file {}?".format(self.frames[0])):
            self.saveas(self.frames[0])
    
    def saveas(self, fpath):
        """Write changed data to molecule/frame file"""
        self.mols.to_file(fpath)
        print("[Info] File {} has been written successfully.".format(fpath))
