import sys
import pyMDMix as pm
iwrap=sys.argv[1]
p = pm.loadProject()
for repl in p.replicas.values():
	print "Setting iwrap %s for relpica %s"%(iwrap, repl.name)
	repl.iwrap=iwrap
	repl.write()

