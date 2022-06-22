import sys
from MDMixTrajectory import Trajectory
import numpy as npy
import Biskit as bi

outprefix=sys.argv[3]
pdb = bi.PDBModel(sys.argv[1])
#pdb = bi.PDBModel('ETA_0/HEWL_ETOH_FREE_ETA_0.pdb')
traj = Trajectory(sys.argv[2],'AMBER',pdb)
#traj = Trajectory('ETA_0/0','AMBER',pdb)
protmask = pdb.maskProtein()
c2mask = pdb.maskFrom('residue_name','ETA')*pdb.maskFrom('name','C2')
c2ids = npy.array(pdb['residue_number'])[c2mask]
cpoint = npy.array([3.850,   2.848 ,  5.845])
tolerance=2.5

print "Parsing trajectory to grab snapshots where ETA is near point %s in %.2f Angstroms"%(cpoint, tolerance)
snapi=0
for nano in traj:
	print "Working on %s"%nano.fname
	for snap in nano:
		snapi += 1
		c2s = snap[c2mask]	
		d=npy.sqrt(((c2s-cpoint)**2).sum(axis=1))
		closeidx=npy.where(d<=tolerance)[0]
		l = npy.take(c2ids,closeidx).tolist()
		if len(l) > 1:
			print "Warning, more than 1 ETA near %i A of cpoint, snap %i"%(tolerance,snapi)
		resmask = pdb.maskFrom('residue_number',l)
		snappdb = pdb.compress(protmask+resmask)
		snappdb.writePdb(outprefix+'_%i.pdb'%snapi)
print "DONE"
