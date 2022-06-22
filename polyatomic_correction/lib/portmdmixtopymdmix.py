#!/usr/bin/env python
# Port old structure to new project format
import sys
import os.path as osp
import pyMDMix
import MDMix.ProjectManager as pm
pyMDMix.setLogger(level='DEBUG')
#pyMDMix.S.SOLVENTDB='/home/daniel/.local/lib/python2.7/site-packages/pyMDMix/data/solventlib/SOLVENTS.db'

if len(sys.argv) < 2:
	sys.exit("Port old mdmix projects to new pyMDMix project format\nUSAGE: python portmdmixtopymdmix.py PATHTOPROJECTFILE\n")

# Give path to MDMix project file
mdmixproj = sys.argv[1]
#mdmixproj='hivpr/hivpr.mproj'

# Load MDMix project
mdmix2 = pm.ProjectManager(mdmixproj)
projpath = mdmix2.project.projectPath
sys.tracebacklimit = 3

# create pyMDMix project
newproj = pyMDMix.Project(name=mdmix2.project.projName+'2')
pyMDMix.browser.goHome()
newproj.createProjectFolder()

# Create pyMDMix.Replicas for each replica
# in mdmix2
restrMask = mdmix2.project.restrainMask
if restrMask == '0': restrMask = 'auto'

ref = osp.join(projpath, mdmix2.project.ref)
for repl in mdmix2:
	print "Porting replica %s"%repl.name

	rpath = osp.join(projpath, repl.replPath)
	path = lambda x: osp.join(rpath, x)
	ss = pyMDMix.SolvatedSystem(name=mdmix2.project.projName, top=path(repl.top), crd=path(repl.crd), pdb=path(repl.pdb), ref=ref)
	sets= pyMDMix.MDSettings(solvent=repl.solvent, restrMask=restrMask, nanos=repl.nanos, temp=repl.temp,
				restrMode=repl.restrMode, restrForce=repl.restrForce,mdoutfiletemplate=repl.productionFileTemplate.replace('nano','step'))
	newrepl = pyMDMix.Replica(system=ss, mdsettings=sets, name=repl.name)
	pyMDMix.browser.goHome()
	pyMDMix.browser.chdir(newproj.replPaths)
	newrepl.createFolder()
	# Import data
	newrepl.importData(minfolder=path(repl.minfolder), eqfolder=path(repl.eqfolder), mdfolder=path(repl.mdfolder),
				alignfolder=path(repl.centerfolder))
				
print "DONE"
