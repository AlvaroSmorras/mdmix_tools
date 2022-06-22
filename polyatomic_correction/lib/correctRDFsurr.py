#! /usr/bin/env python
# Run correction by solvent name or by list of replica names
import os
import os.path as osp
import logging
import warnings
import tempfile
import numpy as npy
import Biskit as bi

#warnings.filterwarnings('error')

import multiprocessing

import pyMDMix
import pyMDMix as pm
import pyMDMix.tools as T
from pyMDMix.Energy import EnergyConversion
from pyMDMix.GridsManager import GridSpace

class CorrectionsError(Exception):
    pass

class CorrectReplicas(object):
    action_name = "CorrectReplicas"
    def __init__(self, replicalist, subregion=None, stepselection=[], outprefix='',
                    protvalue=999., outpath='CORRECTED', ncpus=1, atomsexcluded={}, 
		    *args, **kwargs):
        """
        """
        self.log = logging.getLogger("CorrectReplicas")
        self.replicalist = replicalist
        #self.log.info("Setting up density grid calculation for all heavy atoms in replica %s"%replica.name)

        for r in self.replicalist:
            if not isinstance(r, pyMDMix.Replica): raise CorrectionsError, "replica argument of wrong type."
	    if not r.isAligned(stepselection): raise CorrectionsError, "Cannot correct non-aligned trajectory"

        if not sum([r.isAligned(stepselection) for r in self.replicalist]) == len(self.replicalist):
            raise CorrectionsError, "Cannot calculate density over non-aligned trajectory"
        
        self.solvents = [r.getSolvent() for r in self.replicalist]
        names = set([s.name for s in self.solvents])
        if  len(names) > 1:
            raise CorrectionsError, "All replicas must have same solvent."
        self.solvent = self.solvents[0]
        self.pdbs = {replica.name: replica.getPDB() for replica in self.replicalist}

        if outprefix: self.outprefix = outprefix
        else: self.outprefix = ''

        # Wrok only on subregion?
        if subregion:
            if len(subregion) != 2:
                raise CorrectionsError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
            if len(subregion[0]) != 3 or len(subregion[1] != 3):
                raise CorrectionsError, "subregion argument should have two xyz coordinates: min and max. E.G.: ((x0,y0,z0),(x1,y1,z1))"
        self.subregion = subregion
        
        # To be set un setup()
        self.probes = []
        self.hainfo = {}
        self.dggrids = {}
        self.protvalue = protvalue
        self.outpath = outpath
        self.ncpus = ncpus
	self.atomsexcluded=atomsexcluded
	self.stepselection = stepselection
        self.setup()
    
    def setup(self):
        self.log.info("Setting up parameters...")
        self.setHAinfo()
        self.prepareGrids()
        self.partialDG = self.__prepareResultsContainers()
        f = self.container.getIndexFunction()
        def fx(x):
            r = f(x)
            if not r: return [-1,-1,-1]
            else: return r
        self.indexFunction = fx
	self.prepareResIter()
        # Set subregion function?
        self.log.info("Ready to calculate")

    def prepareResIter(self):
	"Prepare indices to identify residues inside each replica pdb"
	self.replResIter = {}
	for i,rep in enumerate(self.replicalist):
	    pdb = rep.getPDB()
	    pdb._SolvatedPDB__prepareSolventResMasks()
	    resnames = [r.name for r in pdb.solvent.residues]
	    self.replResIter[i]= {r:[npy.where(m)[0] for m in pdb.resMasks[r]] for r in resnames}

    def setHAinfo(self):
        "Prepare a list with all Heavy atoms to track and center of mass for each residue in the solvent"
        for res in self.solvent.residues:
            if res == 'WAT' or res == 'HOH': continue # Don't track water
            excluded = self.atomsexcluded.get(res.name) or []  # Empty list if nothing to exclude for this residue name or return list of atom names to exclude
            if excluded: 
		self.log.info("Excluding atoms %s for residue %s from correction calculation"%(excluded, res))
            self.hainfo[res.name] = dict((a.name,a.id) for a in res.atoms if (a.element != 1 and a.name not in excluded)) # take non-hydrogen atoms only (id and name)
            self.probes.extend(['%s_%s'%(res.name, a.name) for a in res.atoms if (a.element != 1 and a.name not in excluded)])
            self.probes.extend(['%s_%s'%(res.name, 'COM')])

                    
    def prepareGrids(self):
        "Fetch density grids for all heavy atoms in solvent and all replicas, add them and convert to free energies"
        self.log.debug("Fetching HA density grids")
        # Check if replica
        self.dgrids = {}
        for probe in self.probes:
            probe_grids = []
            for r in self.replicalist:
                r.go()
                d = r.getGridsByProbe(probe)
                if d:
                    for g in d[probe]:
                        if g.type == 'MDMIX_DENS': probe_grids.extend([g])
            if len(probe_grids) != len(self.replicalist):
                raise CorrectionsError, "Not all replicas have probe %s density grids"%probe
            self.dgrids[probe] = probe_grids

        # Save one grid as container
        self.container = self.dgrids.values()[0][0].copy()
        self.container.data *= 0

        # Convert density grids to energy
        # Sum densities, calculate expected num for each replica and convert
        self.log.info("Converting density grids to energy...")
        ene = EnergyConversion()
        sumexpectval = {}
	standardstate = {}
        for res in self.solvent.residues:
            if res == 'WAT' or res == 'HOH': continue
            comp = res.name+'_COM'
            expected = [ene.calcReplicaExpectedValue(r, comp) for r in self.replicalist]
            sumexpectval[res.name] = npy.sum(expected)
	    standardstate[res.name] = ene.calcDG0correction(self.replicalist, res.name)

        self.sumdensgrids = {}
        for probe in self.probes:
            self.log.info("Converting %s"%probe)
     	    res = probe.split('_')[0]
            expect = sumexpectval[res]
	    stcorrect = standardstate[res]
            # Averaging if more than 1 replica
            if len(self.replicalist) == 1:
                sumDensityGrid = self.dgrids[probe][0]   # Not doing average!
                self.log.warn("Only one replica for conversion! Can not do average of \
                        the grids. Saving same grid with avg output name. Replica: %s"%self.replicalist[0].name)
                t = self.replicalist[0].temp
            else:
                gspace = GridSpace()
                t = npy.mean([repl.temp for repl in self.replicalist]) # Average temperature for energy conversion
#                    gspace.setT(t)
                gspace.loadGrids(self.dgrids[probe])
                sumDensityGrid = gspace.sum()

            self.sumdensgrids[probe] = sumDensityGrid.copy()

            self.log.info("Expected value: %.3f, Temperature: %.2f"%(expect, t))
            # Convert sumDensityGrid to free energies and apply penalty
            self.RT = t* 0.001987

            # Get a mask of zero values (probably correspoding to protein positions
            maskzeros = sumDensityGrid.data == 0
            sumDensityGrid.data[maskzeros] = 1 # avoid zero-divisions
            dgData = -self.RT*npy.log(sumDensityGrid.data / float(expect))
	    if '_COM' in probe: 
		dgData += stcorrect  # correct by standard state if center of mass
		self.log.info("Correcting COM probe %s by standard state (value: %.3f kcal/mol)"%(probe, stcorrect))

            dgData[maskzeros] = self.protvalue # set back value to identify zeros

            suffix = '_DG'
            sumDensityGrid.update(dgData)
            # Saving to output directory
            if pyMDMix.browser.getcwd() != pyMDMix.browser.home: pyMDMix.browser.goHome()
            if not osp.exists(self.outpath): os.mkdir(self.outpath)
            outname = osp.abspath(self.outpath+os.sep+self.outprefix+probe+suffix+'.dx')
            sumDensityGrid.setType('MDMIX_RAW_AVG')
            sumDensityGrid.setProbe(probe)
            self.log.info("Saving grid %s..."%outname)
            sumDensityGrid.writeDX(outname)
            self.dggrids[probe] = sumDensityGrid

	    # Save cumulative density grids
	    self.sumdensgrids[probe].setType('MDMIX_DENS')
	    self.sumdensgrids[probe].setProbe(probe)
	    outname = osp.abspath(self.outpath+os.sep+self.outprefix+probe+'_cumdens.dx')
	    self.sumdensgrids[probe].writeDX(outname)
	    self.log.info("Saving grid %s..."%outname)

    def __prepareResultsContainers(self):
        """
        Create 1 numpy memmap for each heavy atom density grid. It will store first the radial distribution function
        at each grid point and finally will calculate the partial contribution to the DG0 of the COM energy grid.
        """
        dic = {}
        for probe, dgg in self.dggrids.iteritems():
            tmp = tempfile.mktemp(prefix='mdmix_mmap_')
            cmmap = npy.memmap(tmp, mode='w+', dtype='float32',shape=dgg.data.shape)
            dic[probe] = cmmap
            self.log.debug("Setting memmap name %s for probe %s"%(tmp, probe))
        return dic

    def run(self):
        # Prepare parallel execution if ncpus > 1
        if self.ncpus > 1:
            queue = multiprocessing.Queue(self.ncpus*10)
            processes = self.prepareParallelWorkers(queue, self.ncpus)
            [p.start() for p in processes]
            
        # Process trajectory for all replicas
        snapi = 0
        for repi,replica in enumerate(self.replicalist):
            self.log.info("+"*20)
            self.log.info("\t\t Replica %s "%replica.name)
            self.log.info("+"*20)
            traj = replica.getTrajectory(stepselection=self.stepselection)
            for fil in traj:
                self.log.info("Working on %s trajectory file: %s"%(replica.name, fil.fname))
                for snap in fil:
                    snapi += 1
		    #print "submitting %i for replica %s"%(snapi,replica.name)
                    if self.ncpus > 1:
			queue.put((snap,snapi,repi))
		    else:
                        self.runFrame((snap, snapi, repi))

        # Terminate queues if parallel
        if self.ncpus > 1:
            [queue.put(None) for _ in range(self.ncpus)]
            [p.join() for p in processes]

        self.log.info("+"*20)
        self.log.info("\t\t DONE ")
        self.log.info("+"*20)

    def runFrame(self, frame):
        # correction in given frame
        # select appropriate residues coordinates inside boundaries (or all if not subRegion)
        snapshot, num, repi = frame
        self.log.debug("+ Repl %d Frame %i"%(repi,num))
        # Add counts to corresponding grid
        for comprobe, res in self.solvent.comprobes.iteritems():
            if res.name == 'WAT' or res.name=='HOH': continue
            for ix in self.replResIter[repi][res.name]:
                rxyz = npy.take(snapshot, ix, axis=0)
                comcoord = rxyz.mean(axis=0)
                comEnergy = self.dggrids[comprobe][comcoord]
                if comEnergy is False or comEnergy == 999.: continue
                atomDGs = []
                atomProbe = []
                for atn, ati in self.hainfo[res.name].iteritems():
                    # Reconstructing probenames by residue
                    probe = res.name+'_'+atn
                    atcoords = rxyz[ati -1]
#                    if self.subregionfx: atcoords = self.subregionfx(atcoords)
		    e = self.dggrids[probe][atcoords]
		    if e == 999. or e is False: e = 0
                    atomDGs.append(e)
                    atomProbe.append((probe,atcoords))

                if npy.any(atomDGs == False): 
		    self.log.debug("Could not find energy value for some atom %s"%atomDGs)
		    continue
                atomDGs = npy.array(atomDGs)
		natoms = atomDGs.size
		dgatomContrib = atomDGs - atomDGs.mean() + (1./natoms)*comEnergy
		atomContrib = npy.exp(dgatomContrib/-self.RT)
                for i, probcoord in enumerate(atomProbe):
		    prob, atcoords = probcoord
                    ene = atomContrib[i]
		    idx = self.indexFunction(atcoords)
                    self.partialDG[prob][tuple(idx)] += ene
                    self.log.debug("RES %s PROBE %s ENE %.2f COM_ENE %.3f"%(res, prob, ene, comEnergy))

    def calcResults(self):
        "Obtain final grids that should be written to disk. Average CMsum over counts and substract rawProbeDG - factor*AVG(CMdg)"
        # Once all trajectory has been parsed though run
        # this command should be executed to finally process results
        self.results = {}
        self.log.info("Processing results...")
        for res, atomsinfo in self.hainfo.iteritems():
            if res == 'WAT': continue
            for probe in self.probes:
                if 'COM' in probe: continue
                # Divide finalDG grids by counts grids for each atom --> NOT NEEDED ANYMORE AS ENERGY WAS DIVIDED PREVIOUSLY
		self.container.update(self.partialDG[probe])
		self.container.writeDX('partialDG_%s.dx'%probe)
		protmask = self.sumdensgrids[probe].data == 0
                data = npy.ma.divide(self.partialDG[probe], self.sumdensgrids[probe].data)
		maskzero = data == 0 
		data[maskzero] = 1  # trick to avoid log 0 error
                data = -self.RT*npy.ma.log(data)
		data[maskzero] = 0.0 # restore zeros
                data.fill_value = self.protvalue # Change NaNs or Inf to zeros
		data = data.filled()
		data[protmask] = self.protvalue
                self.container.update(data)
		# remove edge
		self.container.contract(0.5)
		self.container.expand(0.5)
#                self.containerGrid.update(self.partialDG[res][atomnames[0]])
                self.container.source = "RDFSurrogateCorrection on residue %s probe %s"%(res, probe)
                self.container.setProbe(probe)
                self.container.setType('MDMIX_CORR')
                self.results[probe] = self.container.copy()
        return self.results

    def prepareParallelWorkers(self, dataQueue, nworkers):
        "Will return as many Processes as probes to study. Each will need a snapshot to run. So when sending snapshots to the task queue, repeat each snap as needed"
        # Substitute grids' data by memmap array/files
        # Instantiate workers with memmap arguments and probe info
        workerList = []
	i=0
        for _ in range(nworkers):
	    i+=1
	    print "Preparing Worker",i
            workerList.append(CorrectWorker(dataQueue, self.solvent.comprobes, self.dggrids, 
					    self.hainfo, self.partialDG, self.replResIter, 
					    self.indexFunction, self.RT))
        return workerList


class CorrectWorker(multiprocessing.Process):
    def __init__(self, dataQueue, comprobes, dggrids, hainfo, partialDG, replResIter, indexFunc, RT):
        multiprocessing.Process.__init__(self)

        self.comprobes = comprobes
        self.dggrids = dggrids
        self.hainfo = hainfo
        self.dataQueue = dataQueue
        self.partialDG = partialDG
        self.replResIter = replResIter
	self.indexFunction = indexFunc
	self.RT = RT

    def run(self):
        while True:
            frame = self.dataQueue.get()
            if frame is None: break  # Poison pill exit processor
            snapshot, num, repi = frame         #second index is frame number 
            for comprobe, res in self.comprobes.iteritems():
                if res.name == 'WAT' or res.name=='HOH': continue
                for ix in self.replResIter[repi][res.name]:
		    rxyz = npy.take(snapshot,ix,axis=0)
                    comcoord = rxyz.mean(axis=0)
                    comEnergy = self.dggrids[comprobe][comcoord]
                    if comEnergy is False or comEnergy == 999.: continue

                    atomDGs = []
                    atomProbe = []
                    for atn, ati in self.hainfo[res.name].iteritems():
                        # Reconstructing probenames by residue
                        probe = res.name+'_'+atn
                        atcoords = rxyz[ati-1]
    #                    if self.subregionfx: atcoords = self.subregionfx(atcoords)
		        e = self.dggrids[probe][atcoords]
		    	if e == 999. or e is False: e = 0
                    	atomDGs.append(e)
                        atomProbe.append((probe,atcoords))

                    if npy.any(atomDGs == False): 
			print "Atom is out: %s %s %.3f %s"%(atomDGs, atomProbe, comEnergy, comcoord)
			continue
                    atomDGs = npy.array(atomDGs)
                    natoms = atomDGs.size
		    dgatomContrib = atomDGs - atomDGs.mean() + (1./natoms)*comEnergy
		    atomContrib = npy.exp(dgatomContrib/-self.RT)
                    for i, probcoord in enumerate(atomProbe):
			prob,atcoords = probcoord
                        ene = atomContrib[i]
			idx = self.indexFunction(atcoords)
			o = self.partialDG[prob][idx]
			self.partialDG[prob][idx] += ene
			#print "Adding prob: %s coord: %s idx: %s old: %.3f ene: %.3f new:%.3f"%(prob, atcoords, idx, o, ene, self.partialDG[prob][idx])
                       

        return

def CorrectReplicas_postprocess(results, outfolder='CORRECTED', prefix='RDFsurr_', **kwargs):
    """
    Save results of a density calculation for replica *replica*.

    :arg dict densityResults: Result from a :class:`Action.Density` call.
    :arg replica: Replica were to save results.
    :type replica: :class:`Replicas.Replica`
    """
    pyMDMix.browser.goHome()
    if not osp.exists(outfolder): os.mkdir(outfolder)
    print results
    for probe, grid in results.iteritems():
        grid.writeDX(osp.join(outfolder, prefix+probe+'.dx'))
    pyMDMix.browser.goback()
    


if __name__=='__main__':
	import sys
	import time

	proj = pm.loadProject()
	pm.setLogger('DEBUG')
	stepselection=range(1,21)

	if len(sys.argv) == 4:
		option = sys.argv[1].lower()
		ncpus = int(sys.argv[2])
		solvent = sys.argv[3].upper()
		replicas = [r for r in proj.replicas.values() if r.solvent == solvent]
	elif len(sys.argv) > 4:
		option = sys.argv[1].lower()
		ncpus = int(sys.argv[2])
		rnames = sys.argv[3:]
		replicas = [r for r in proj.replicas.values() if r.name in rnames]
	else:
		sys.exit("USAGE: 1) correctRDF.py [DENSITY|CORRECT|ALL] NCPUS SOLVENT\n 2) correctRDF.py [DENSITY|CORRECT|ALL] NCPUS REPLNAME REPLNAME ...")

	t0 = time.time()

	if option == 'density' or option == 'all':
		print "CALCULATING DENSITY FOR ALL HA"
		actions = pm.Analysis.ActionsManager(ncpus)
		actions.addReplicas(replicas)
		actions.addActions('DensityGridsAllHA')
		actions.prepareRun(stepselection=stepselection)
		actions.run(stepselection=stepselection)
		actions.processResults()

	elif option == 'correct' or option == 'all':
		print "CONVERTING TO ENERGIES AND CORRECTING"
		correct = CorrectReplicas(replicas, ncpus=ncpus, atomsexcluded={'MAM':['C2'],'COO':['C1'],'TRI':['F1','F2','F3'],'MSU':['S1'],'CLE':['LP1']}, stepselection=stepselection)
		correct.run()
		results = correct.calcResults()
		CorrectReplicas_postprocess(results)

	else:
		print "INVALID OPTION. SHOULD BE: DENSITY or CORRECT or ALL"

	print 'DONE in ', t0-time.time()
