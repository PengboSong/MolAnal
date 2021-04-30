# coding=utf-8

import os


class IOMatrix(object):
    # PDB template line
    PDB_FORMAT = ("HETATM{0:>5} {1:<4} {2:>3}  {3:>4}    {4:>8.3f}{5:>8.3f}{6:>8.3f}"
                  "{7:>6.2f}{8:>6.2f}          {9:>3}\n")
    # GRO template line
    GRO_FORMAT = "{0:>5}{1:<4}  {2:>4}{3:>5}{4:>8.3f}{5:>8.3f}{6:>8.3f}"
    GRO_VELOCITY = "{0:>8.3f}{1:>8.3f}{2:>8.3f}"
    # XYZ template line
    XYZ_FORMAT = "{0:>8.3f}{1:>8.3f}{2:>8.3f}"

    def from_pdb(self, fpath):
        """Read data from PDB format file"""
        raise NotImplementedError('Reading from PDB format file not supported.')

    def from_gro(self, fpath):
        """Read data from GRO format file"""
        raise NotImplementedError('Reading from GRO format file not supported.')

    def from_xyz(self, fpath):
        """Read data from XYZ format file"""
        raise NotImplementedError('Reading from XYZ format file not supported.')

    def from_file(self, fpath):
        """Read data from file"""
        suffix = os.path.splitext(fpath)[-1]
        if os.path.isfile(fpath):
            if suffix == '.pdb':
                self.from_pdb(fpath)
            elif suffix == '.gro':
                self.from_gro(fpath)
            elif suffix == '.xyz':
                self.from_xyz(fpath)

    def to_pdb(self, fpath):
        """Write data to PDB format file"""
        raise NotImplementedError('Writing to PDB format file not supported.')

    def to_gro(self, fpath):
        """Write data to GRO format file"""
        raise NotImplementedError('Writing to GRO format file not supported.')

    def to_xyz(self, fpath):
        """Write data to XYZ format file"""
        raise NotImplementedError('Writing to XYZ format file not supported.')

    def to_file(self, fpath):
        """Write data to file"""
        suffix = os.path.splitext(fpath)[-1]
        if suffix == '.pdb':
            self.to_pdb(fpath)
        elif suffix == '.gro':
            self.to_gro(fpath)
        elif suffix == '.xyz':
            self.to_xyz(fpath)
