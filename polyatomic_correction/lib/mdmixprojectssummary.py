import sys
import os
import fnmatch
import pyMDMix as pm

availprojects = []
cwd = os.getcwd()

for root, dirs, files in os.walk('.'):
	projs = fnmatch.filter(files, '*mproj')
	if projs:
		for p in projs:
			ppath = os.path.join(root,p)
			try:
				pm.browser.chdir(cwd)
				print cwd
				availprojects.append(pm.loadProject(ppath))
				print "Found project file ",ppath
			except:
				print "Old MDMix project found: ",ppath				
	
print 'done'
