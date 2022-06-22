# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$01-jun-2011 15:31:10$"

import sys
import os
import glob

import Biskit as bi
from pyMDMix.NamdDCDParser import NamdDCDParser

class TrajFile:
    def __init__(self, file, pdb, prog):
        "pdb is a PDBModel, prog is either 'amber' or 'namd'."
        self.fname = file
        self.pdb = pdb
        self.traj = None
        self.prog = prog
        self.nframes = 0
        self.loadFile()

    def loadFile(self):
        if self.prog == 'amber':
            self.traj = bi.AmberCrdParser(self.fname, self.pdb.fileName, 1)
            self.traj.crd.readline()    #skip first line
        else:
            #NAMD
            from pyMDMix.NamdDCDParser import NamdDCDParser
            self.traj = NamdDCDParser(self.fname, self.pdb, box=1)

    def nextFrame(self):
        try:
            if self.prog == 'amber':
		frame = self.traj.nextFrame()
            else: frame = self.traj.read_dcdstep()
        except:
	    print 'frame did not load properly or was last frame (%s)' %self.nframes
            return False

        self.nframes += 1
        return frame
    def close(self):
        if self.prog == 'namd': self.traj.close()
        else: pass

class Trajectory:
    """Allow reading trajectories from NAMD or AMBER giving a directory
        with the files in it as 'path' argument."""
    def __init__(self, path, mdprog, pdb, ext=False):
        self.path = path
        self.mdprog = mdprog.lower()
        self.ext = ext
        if isinstance(pdb, str): self.pdb = bi.PDBModel(pdb)
        elif isinstance(pdb, bi.PDBModel): self.pdb = pdb
        else: print "WARNING! Unkown pdb argument type."
        self.files = self.getFileList()
        self.nfiles = len(self.files)
    def getFileList(self):
        if self.ext:
            ext = self.ext
        elif self.mdprog == 'amber':
            ext = 'x'
        elif self.mdprog == 'namd':
            ext = 'dcd'
        else:
            sys.exit("ERROR with program name. Should be either namd or amber")

        print "Will use .%s extension for trajectory files."%ext
        flist = glob.glob(self.path+os.sep+'*.'+ext)
	print(self.path+os.sep+'*.'+ext)
        if len(flist) == 0:
            sys.exit("No files with extension .%s found in path %s"%(ext,self.path))

        return flist

    def readFile(self, i):
        "i is the index of the file in the files list"
        self.file = TrajFile(self.files[i], self.pdb, self.mdprog)

if __name__ == "__main__":
    print "Hello! This is not an executable module"
