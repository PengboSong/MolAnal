# coding=utf-8

import math
import os

import numpy as np
from custom.general import listFiles

from gmx.math.common import fitplane
from gmx.other.command_line import ConsoleBase
from gmx.other.data_type import GMXDataType
from gmx.other.mol_order import MolOrder
from gmx.other.register_func import RegisterFunction
from gmx.structure.mol_matrix import MolMatrix


class MolEdit(ConsoleBase, RegisterFunction):
    PROMPT = "Molecules Editing for GROMACS"
    EXIT_PROMPT = "Press any key to exit..."
    PLACEHOLDER = "[MolEdit]>>> "

    def __init__(self, frame=''):
        super().__init__()
        self.molmat = MolMatrix()
        self.frames = []
        if os.path.isfile(frame):
            self.molmat.from_file(frame)
            self.frames.append(frame)

        # Register functions
        self.register_("load_mol", self.load_mol, str)
        self.register_("copy_mol", self.copy_mol, int)
        self.register_("move", self.move, int, list)
        self.register_("moveto", self.moveto, int, list)
        self.register_("orrient", self.orrient, int, list)
        self.register_("random_orrient", self.random_orrient, int)
        self.register_("listmols", self.listmols)
        self.register_("resort", self.resort)
        self.register_("save", self.save)
        self.register_("saveas", self.saveas, str)
        self.register_("help", self.help_)
        self.register_("exit", self.exit_)

        self.register_("genbox", self.molmat.genbox, [int, float, list])
        self.register_("extend-z", self.molmat.extend_z, [int, float])
    
    def __call__(self):
        """Initialize command-line interactive interface"""
        print(self.PROMPT)
        # Read commands from console
        while not self.exit_signal_:
            self.launch_(input(self.PLACEHOLDER))
            #try:
            #    self.launch_(input(self.PLACEHOLDER))
            #except Exception as err:
            #    print("[Error] {}".format(err))
        if len(self.frames) > 0:
            if self.askyn(prompt="Save modifications?"):
                self.save()
    
    def check_molid(self, molid):
        if molid > 0 and molid <= len(self.molmat):
            return True
        else:
            return False

    def load_mol(self, frame=''):
        """Load another molecule/frame file into current molecule matrix.
        
        Args:
            frame(str): Path to molecule/frame file.
        """
        if os.path.isfile(frame):
            new_mol = MolMatrix()
            new_mol.from_file(frame)
            self.molmat.merge(new_mol)
            self.frames.append(frame)
    
    def copy_mol(self, copymol):
        """Copy a selected molecule.
        
        Args:
            copymol(int): Index of molecule to copy (start from 1).

        Raise:
            ValueError: Molecule not found.
        """
        # Check
        if not self.check_molid(copymol):
            raise ValueError("Target molecule %d to cut not found." % copymol)
        newmol = self.molmat[copymol - 1].copy()
        self.molmat.append(newmol)
    
    def cut(self, centermol, shape="sphere", parameter=(1.0,)):
        """Cut molecules near central molecule within a custom shape."""
        # Check
        if not self.check_molid(centermol):
            raise ValueError("Target molecule %d to cut not found." % centermol)
        # Calculate coordinate of target molecular center
        center = self.molmat[centermol - 1].center()
        # Cut a sphere around target molecular center
        if shape == "sphere":
            inside = MolMatrix()
            outside = MolMatrix()
            radius = parameter[0]
            # Scan all molecules
            for mol in self.molmat:
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
        self.molmat[movemol - 1].move(vector)

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
        movevec = point - self.molmat[movemol - 1].center()
        self.molmat[movemol - 1].move(movevec)
    
    def orrient(self, rotmol, vector):
        """Rotate molecule to a specified direction.
        
        Args:
            rotmol(int): Index of molecule to rotate (start from 1).
            axis(list[int*3]): Vector [vx, vy, vz] as the rotation direction.

        Raise:
            ValueError: Molecule not found.
        """
        # Check
        vector = np.asarray(vector, dtype=GMXDataType.REAL).reshape(-1)
        assert vector.size == 3, "Direction vector should be like [ax, ay, az]."
        if not self.check_molid(rotmol):
            raise ValueError("Target molecule %d to rotate not found." % rotmol)
        # Normalize
        vector /= np.linalg.norm(vector)
        # Get the normal vector determined by molecular coordinates
        plparms = fitplane(self.molmat[rotmol - 1].xyzs)[:3]
        normal = np.asarray(plparms, dtype=GMXDataType.REAL)
        normal /= np.linalg.norm(normal)
        # Rotation angle
        cosang = min(1., max(-1., np.dot(normal, vector)))   # Cosine in range [-1., 1.]
        sinang = math.sqrt(1 - cosang**2)
        if (1. - cosang) < 1e-3:   # No need to rotate
            pass
        else:
            # Rotation axis
            axis = np.cross(normal, vector)
            # Rotate
            self.molmat[rotmol - 1].rotate(axis, cosang, sinang)

    def random_orrient(self, rotmol):
        """Rotate molecule to a random direction.
        
        Args:
            rotmol(int): Index of molecule to rotate (start from 1).

        Raise:
            ValueError: Molecule not found.
        """
        self.orrient(rotmol=rotmol, axis=np.random.random(3, dtype=GMXDataType.REAL))
    
    def listmols(self):
        """List all molecules in system"""
        for molnm, molids in self.molmat.clustering().items():
            print('- {} : {}'.format(molnm, ', '.join(str(i) for i in molids)))

    def resort(self):
        """Resort molecules in molecular name order"""
        self.molmat.resort(MolOrder.MOL_ORDER)
    
    def save(self):
        """Save to the first loaded molecule/frame file"""
        if self.askyn(prompt="Confirm overwrite file {}?".format(self.frames[0])):
            self.saveas(self.frames[0])
    
    def saveas(self, fpath):
        """Write changed data to molecule/frame file"""
        self.molmat.to_file(fpath)
        print("[Info] File {} has been written successfully.".format(fpath))
