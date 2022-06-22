# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$01-jun-2011 15:14:16$"

import sys
import os
import time
import copy
import ConfigParser

import numpy as npy
import Biskit as bi

import pyMDMix.GridData as gr

from parseTrajectory import Trajectory


def getOptions(filename):
    config = ConfigParser.ConfigParser()
    config.read(filename)

    # Solvent should be either MAM or ETA
    fileSection = dict(config.items('OPTIONS'))
    solvent = fileSection.get('solvent').upper()
    if solvent != 'MAM' and solvent != 'ETA' and solvent != 'ION':
        print "ERROR in configuration file. SOLVENT should be ETA, MAM or ION."
        sys.exit(1)
        
    # MD program must be either NAMD or AMBER
    mdprog = fileSection.get('mdprog').upper()
    if mdprog != 'NAMD' and mdprog != 'AMBER':
        print "ERROR in configuration file. MDPROG should be NAMD or AMBER"
        sys.exit(1)
    
    # Outpath
    outpath = fileSection.get('outpath')
    if not outpath:
        outpath = './'
        
    # Trajectory path
    trajpath = fileSection.get('trajpath')

    # PDB of the solvated system
    syspdb = fileSection.get('syspdb')

    # Trajectory file extension
    ext = fileSection.get('ext')
    if ext:
        print "Will use trajectory file extension ."+ext

    # DG grids
    fileSection = dict(config.items('GRIDS'))
    grids_dict = dict(zip([k.upper() for k in fileSection.keys()],fileSection.values()))

    opts = [solvent, mdprog, outpath, trajpath, syspdb, grids_dict, ext]
    return opts

class CorrectGrid:
    
    def __init__(self, type, atom, raw_grid, ASAfactor=1):
        self.type = type
        self.atom = atom
        self.raw = raw_grid
        self.asa = ASAfactor
        self.counts = npy.zeros_like(self.raw.data) if type != 'CM' else None
        self.corrected = npy.zeros_like(self.raw.data) if type != 'CM' else None
        
    def correctPoint(self, XYZ, CMDG):
        point = self.raw.getIndex(XYZ)
        #print(point, self.raw.data[point])
        if point and self.raw.data[point] < 999.:
            self.counts[point] += 1
            #print(xyz, self.raw.data[point], self.asa*CMDG)
            self.corrected[point] += self.raw.data[point] - self.asa*CMDG
        else:
            print "Point out of grid:",self.atom,XYZ
    def averageByCount(self):
        self.corrected /= self.counts
        self.corrected[npy.isnan(self.corrected)] = 0.0

    def _write(self, data, name=False, prefix=''):
        if name:
            outname = name
        else:
            sourcefile = self.raw.source
            outname = prefix+os.path.basename(sourcefile)
        g = copy.deepcopy(self.raw)
        g.data = data
        print "Writing Grid: ",outname
        g.writeDX(outname)
        #npy.save(outname+'.npy', data)

    def writeCorrectedGrid(self, name=False, prefix='corrected_'):
        self._write(self.corrected, name=name, prefix=prefix)

    def writeCountGrid(self, name=False, prefix='counts_'):
        self._write(self.counts, name=name, prefix=prefix)


class CorrectMolecule:
    valid_solvents = ['ETA','MAM','COO','CN3']
    solvent_atoms = {'ETA':{'CT':'C1','OH':'O1','CM':'C2'},
                     'MAM':{'O':'O1','N':'N1','CM':'C2'},
                     'CN3':{'N':'N1', 'CM_N':'C1'},
                     'COO':{'O':['O1','O2'],'CM_O':'C2'}}
    ASA_factors = {'ETA':{'CT':0.55648,'OH':0.72549,'CM':1},
                    'MAM':{'N':0.67891,'O':0.777314,'CM':1},
                    'COO':{'O':0.48719,'CM':1},
                    'CN3':{'N':0.55650,'CM':1}}
                    
    def __init__(self, solvent):
        self.solvent = solvent.upper() if solvent.upper() in CorrectMolecule.valid_solvents else sys.exit('Invalid solvent')
        self._prepare()
    def _prepare(self):
        # Prepare atom names and types according to the solvent name
        self.types_atoms = CorrectMolecule.solvent_atoms[self.solvent]

        # Prepare atom names
        self.atoms = []
        for t in self.types_atoms.values():
            if isinstance(t,list):
                [self.atoms.append(el) for el in t]
            else:
                self.atoms.append(t)
                
        # Types
        self.types = self.types_atoms.keys()

        # Map names to types
        # check if one type has multiple atom names
        a_to_t = []
        for t,n in self.types_atoms.iteritems():
            if t == 'CM_N' or t == 'CM_O': t='CM'
            if isinstance(n, list):
                [a_to_t.append([el,t]) for el in n]
            else:
                a_to_t.append([n,t])
        self.atoms_types = dict(a_to_t)
        
        # If ION, change CM_N and CM_O to CM
#        if self.solvent == 'COO':
#            self.types_atoms['CM'] = self.types_atoms['CM_O']
#            self.types.pop(self.types.index('CM_O'))
#            self.types.append('CM')
#        if self.solvent == 'CN3':
#            self.types_atoms['CM'] = self.types_atoms['CM_N']
#            self.types.pop(self.types.index('CM_N'))
#            self.types.append('CM')
            
        # Others
        self.grids = {}
        self.averaged = False
        
    def addGrids(self, grids):
        "Parse raw DG grids as a dictionary using the corresponding grid atomtype\
        as key, and the grid instance as value. From the dict, only needed grids will be added"
        # Check the dictionary contains the grids expected
        filtered_grids = dict([(t,g) for t,g in grids.iteritems() if t in self.types])
        
        if len(filtered_grids.keys()) != len(self.types):
            allgrids = [g not in grids.keys() for g in self.types]
            miss = npy.array(self.types)[allgrids]
            print "addGrids ERROR: Missing grids:",miss
            sys.exit(1)
            
        # Check all values are Grid instances
        is_grid = [not isinstance(g, gr.GridData) for g in filtered_grids.values()]
        if npy.any(is_grid):
            print "addGrids ERROR: grids dictionary parameter should contain Grid instances as values"
            
        # If all is correct assign each grid to each key and atom name
        # CM grid will be stored separately
        for type, grid in filtered_grids.iteritems():
                if type == 'CM_O' or type == 'CM_N': newtype = 'CM'
                else: newtype = type

                # If grid is for CM, set all values > 0 to zero
                if newtype == 'CM':
                    grid.data[grid.data > 0] = 0.0
                    
                self.grids[newtype] = CorrectGrid(newtype, self.types_atoms[type], grid,
                ASAfactor=CorrectMolecule.ASA_factors[self.solvent][newtype])
#            else:
#                self.CM = CorrectGrid(type, self.types_atoms[type], grid)

    def correctMolecule(self, xyz, atoms):
        "xyz is a npy.array with all the coordinates of a snapshot for a given molecule\
        atoms is a list with the atom names to identify the xyz array coordinates"
        if len(xyz) != len(atoms):
            print "Invalid number of types or xyz coordinates"
        types = [self.atoms_types[at] for at in atoms]
#        print xyz, atoms, types
        coords = dict(zip(types,xyz))
        CMgrid = self.grids['CM'].raw
        CMpoint = CMgrid.getIndex(coords['CM'].tolist())
	#print(CMpoint,CMgrid.data[CMpoint])
        if CMpoint and CMgrid.data[CMpoint] < 999.:
            CMDG = CMgrid.data[CMpoint]
            for t,p in coords.iteritems():
                if t != 'CM':
                    self.grids[t].correctPoint(p, CMDG)
        else:
#            print "CM out of grid!"
            pass
            
    def averageGrids(self):
        if not self.averaged:
            for t,g in self.grids.iteritems():
                if t != 'CM':
                    self.grids[t].averageByCount()
            self.averaged = True

    def writeCorrectedGrids(self, outpath="./", prefix='corrected_'):
        if not self.averaged:
            self.averageGrids()
        
        for t,g in self.grids.iteritems():
            if t != 'CM':
                self.grids[t].writeCorrectedGrid(prefix=outpath+prefix)

    def writeCountGrids(self, outpath='./', prefix='counts_'):
        for t,g in self.grids.iteritems():
            if t != 'CM':
                self.grids[t].writeCountGrid(prefix=outpath+prefix)
                
if __name__ == "__main__":

    configfile = sys.argv[1]
    solvent, mdprog, outpath, trajpath, syspdb, grids_dict, extension = getOptions(configfile)
    
    # Prepare the grid files
    # and check they all have same shape and origin
    print "Loading DG grids"
    shapes = []
    origins = []
    for type, grid in grids_dict.iteritems():
        grids_dict[type] = gr.get_grid(grid)
        shapes.append(grids_dict[type].data.shape)
        origins.append(tuple(grids_dict[type].origin))
        
    if len(npy.unique(shapes)) != 1:
        sys.exit("ERROR: Shape missmatch between grids.")
    if len(npy.unique(origins)) != 1:
        sys.exit("ERROR: Origin missmatch between grids.")
            
    # Create instance/s of CorrectMolecule
    # if ION simulation, two instances must be created (COO and CN3)
    correctSolvList = []
    if solvent == 'ION':
        correctSolvList.append(CorrectMolecule('COO'))
        correctSolvList.append(CorrectMolecule('CN3'))
    else:
        correctSolvList.append(CorrectMolecule(solvent))
    
    # Add grids.
    [s.addGrids(grids_dict) for s in correctSolvList]
    # Load the PDBfile and prepare mask
    # for the atoms we are interested in
    print "Loading PDB and preparing masks and residue identification"
    pdb = bi.PDBModel(syspdb)
    #pdb.fixNumbering(force=True)
    resnums = npy.array(pdb['residue_number'])

    solv_masks_list = []
    for s in correctSolvList:
        solv_masks_list.append(pdb.maskFrom('residue_name',s.solvent) * pdb.maskFrom('name',s.atoms))
        
    # Start parsing trajectory
    AllTraj = Trajectory(trajpath, mdprog, pdb, ext=extension)
    for traj_i in range(AllTraj.nfiles):
        print "Working on file ",AllTraj.files[traj_i]
        AllTraj.readFile(traj_i)
        print(traj_i)
        # Loop over frames in file
        frame = AllTraj.file.nextFrame()
        fi = 1
        while isinstance(frame, npy.ndarray):
            print "\t+ Snapshot ",fi
            
            # Loop over different solvent instances to correct and their masks
            for solv, solv_mask in zip(correctSolvList,solv_masks_list):

                solvframe = frame[solv_mask,]
                resids = resnums[solv_mask]
                u_resids = npy.unique(resids)
                atnames = npy.array(pdb['name'])[solv_mask]

                for res in u_resids:
                    ids = npy.where(resids == res)[0]
                    xyz = npy.take(solvframe, ids, axis=0)
                    names = npy.take(atnames, ids)
                    solv.correctMolecule(xyz, names)
                
            fi += 1
            frame = AllTraj.file.nextFrame()
    
    for cSolv in correctSolvList:
        cSolv.writeCorrectedGrids(outpath=outpath, prefix='poly_correctedDG_')
#        cSolv.writeCountGrids(outpath=outpath, prefix='counts_')
