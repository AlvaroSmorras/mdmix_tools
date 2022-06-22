#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="dalvarez"
__date__ ="$04-abr-2011 16:53:13$"

import math,sys
import numpy as npy
import Biskit as bi
import pygraph
from pygraph.classes.graph import graph
#from pygraph.algorithms.heuristics import euclidean
#from pygraph.algorithms.heuristics import chow
import pygraph.algorithms.minmax as minmax
#from pygraph.readwrite.dot import write
#import gv
AdditionError = pygraph.classes.exceptions.AdditionError
#import Biskit.BiskitFF as bff

try:
    import Biskit.cudm as cuda
    _cudaAvail = True
except ImportError:
    _cudaAvail = False

def pointsOnSphere(N):
    N = float(N) # in case we got an int which we surely got
    pts = []

    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / N
    for k in range(0, int(N)):
        y = k * off - 1 + (off / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        pts.append([math.cos(phi)*r, y, math.sin(phi)*r])

    return npy.array(pts)

def pos_distance_sq(point1, point2):
  return npy.sum((point2-point1)**2)


def buildVdw(pdb):
    "generates a BiskitProfile with he VdW radii for each atom extracted from Amber topology"
    rdict = {}
    vdwr = []
    for at in pdb:
        '''
        type = at['amber_atype']
        if not rdict.has_key(type):
            r = bff.getEquilibriumDistance(pdb, type, type, verbose=True)
            d = r['distance']/2.
            if not npy.isnan(d): rdict[type] = d
            else: rdict[type] = 0
        vdwr.append(rdict[type])
        '''
        vdwr.append(at['alternate'])
    pdb['vdwradii'] = npy.array(vdwr)


def calculateASAamber(pdb,  probe, n_sphere_point=960, getpoints=False, radii=False, include_list=None):
    """
    Calculate atomic ASA. pdb is a PDBModel with the Amber Topology loaded.
    getpoints True if you want the 3D point coordinates on the surface calculated
    radii True if you want for each ASA, the radius of the atom (+ probe) aswell
    If include_list with indices is given, only that atoms are taken into the calculation
    """
    buildVdw(pdb)
    vdwradii = pdb['vdwradii'] + probe
    sphere_points = pointsOnSphere(n_sphere_point)
    distMat = calculateDM(pdb.xyz, pdb.xyz)
    nhood = distMat < 8*8
    const = 4.0 * math.pi / len(sphere_points)
    if getpoints: res = []
    areas = []
    totalats = len(pdb)
    for i, r in enumerate(vdwradii):
#        sys.stdout.write("Calculating ASA for atom %i of %i \r" % (i+1,totalats))
        atom_pos = pdb.xyz[i]
        neigh_mask = nhood[i,:]
        neigh_mask[i] = False
        neigh_xyz = pdb.xyz[neigh_mask]
        neigh_r = vdwradii[neigh_mask]
        test_points = sphere_points*r + atom_pos
        inside_neigh = calculateDM(neigh_xyz, test_points) - (neigh_r*neigh_r)
        accessible_m = npy.all(inside_neigh > 0,axis=1)
        num_accessible = accessible_m.sum()

        if getpoints and num_accessible > 0:
            [res.append(point.tolist()) for point in test_points[accessible_m]]
        area = const*num_accessible*r*r
        if radii: areas.append([area, r])
        else: areas.append(area)

    if getpoints:
        return areas, res
    else:
        return areas

def calculateASAamber2(xyz, vdw, probe, n_sphere_point=960, getpoints=False, radii=False, include_list=None):
    """
    Calculate atomic ASA.

    xyz are the protein coordinates
    vdw is a numpy array with the vdw radii for each atom in the protein (same order as xyz)
    
    getpoints True if you want the 3D point coordinates on the surface calculated
    radii True if you want for each ASA, the radius of the atom (+ probe) aswell
    If include_list with indices is given, only that atoms are taken into the calculation
    """
    
    vdwradii = vdw + probe
    sphere_points = pointsOnSphere(n_sphere_point)
    distMat = calculateDM(xyz, xyz)
    nhood = distMat < 8*8
    const = 4.0 * math.pi / len(sphere_points)
    if getpoints: res = []
    areas = []
    totalats = len(vdw)
    for i, r in enumerate(vdwradii):
#        sys.stdout.write("Calculating ASA for atom %i of %i \r" % (i+1,totalats))
        atom_pos = xyz[i]
        neigh_mask = nhood[i,:]
        neigh_mask[i] = False
        neigh_xyz = xyz[neigh_mask]
        neigh_r = vdwradii[neigh_mask]
        test_points = sphere_points*r + atom_pos
        inside_neigh = calculateDM(neigh_xyz, test_points) - (neigh_r*neigh_r)
        accessible_m = npy.all(inside_neigh > 0,axis=1)
        num_accessible = accessible_m.sum()

        if getpoints and num_accessible > 0:
            [res.append(point.tolist()) for point in test_points[accessible_m]]
        area = const*num_accessible*r*r
        if radii: areas.append([area, r])
        else: areas.append(area)

    if getpoints:
        return areas, res
    else:
        return areas


def printXYZ(xyz, name):
    o = open(name, 'w')
    o.write(str(len(xyz))+"\n\n")
    [o.write("C %f %f %f\n"%(r[0],r[1],r[2])) for r in xyz]
    o.write("\n")
    o.close()

def calcContactVolume(area1, area2):
    """
    Calculate the spherical cone volume between 2 areas of two spheres with
    different radii

    see http://mathforum.org/library/drmath/view/55253.html for sphere forumlas

    volume of a truncated cone: V = 1/3*pi*(r1^ + r1*r2 + r2^2)*h
    Thus if r2 is zero (no surface) its the formula of a normal cone.

    """
    # check A1 radii is smaller than A2
    # or swap them
    if area1[1] < area2[1]:
        A1, R1 = area1
        A2, R2 = area2
    else:
        A1, R1 = area2
        A2, R2 = area1

    const = (1/3.)*math.pi

    # small radii cap volume and slice radii (c1)
    h1 = A1/(2*math.pi*R1)
    h1sq = h1**2
    rh = R1*h1
    hcone1 = R1 - h1
    c1 = npy.sqrt((2*rh) - (h1sq))
    Vcap1 = const*h1*(3*rh - h1sq)

    # big radii cap volume, and slice radii (c2)
    h2 = A2/(2*math.pi*R2)
    h2sq = h2**2
    rh = R2*h2
    hcone2 = R2 - h2
    c2 = npy.sqrt((2*rh) - (h2sq))
    Vcap2 = const*h2*(3*rh - h2sq)

    # cone volume (c1, c2, and h of the cone is needed)
    hcone = hcone2 - hcone1
    Vcone = const*(c1**2 + c1*c2 + c2**2)*hcone

    # reconstruct the volumes by: Vcone - Vcap1 + Vcap2
    Vtotal = Vcone - Vcap1 + Vcap2
    
    return Vtotal

def getContactVolumes(PDBModel, probe1, probe2, points1=150, points2=250):
    "Calculates the volume of the contact space for each atom on a biskit.PDBModel\
    probe1 and probe2 are the probe radii for ASA calculation (floats), points1 and points2\
    are related to the accuracy of the ASA calculation (number of sphere points)"
    
    print "Calculating atomic ASA for probe ",probe1
    a1 = calculateASAamber(PDBModel, probe1, n_sphere_point=points1, radii=True)
    print "Calculating atomic ASA for probe ",probe2
    a2 = calculateASAamber(PDBModel, probe2, n_sphere_point=points2, radii=True)
    print "Calculating volumes"
    volumes = [calcContactVolume(a1[i], a2[i]) for i in range(len(PDBModel))]
    serial = PDBModel['serial_number']

    return dict(zip(serial,volumes))

def getContactVolumes2(xyz, vdw, serialn, probe1, probe2, points1=150, points2=250):
    "Calculates the volume of the contact space for each atom on a biskit.PDBModel\
    probe1 and probe2 are the probe radii for ASA calculation (floats), points1 and points2\
    are related to the accuracy of the ASA calculation (number of sphere points)"

    if probe1 < 1.4:
        print "Small probe ASA calculation. Identifying core atoms to prevent intersticial areas"
        a0 = calculateASAamber2(xyz,vdw,1,n_sphere_point=200)
        mask_core = npy.array(a0) < 0.005
        vdw_small = vdw[:]
        vdw_small[mask_core] += 0.5
    else:
        vdw_small = vdw
    print "Calculating atomic ASA for probe ",probe1
    a1, xyz1 = calculateASAamber2(xyz, vdw_small, probe1, n_sphere_point=points1, getpoints=True, radii=True)
    print "Calculating atomic ASA for probe ",probe2
    a2, xyz2 = calculateASAamber2(xyz, vdw, probe2, n_sphere_point=points2, getpoints=True, radii=True)
    print "Calculating volumes"
    volumes = [[calcContactVolume(a1[i], a2[i]),a1[i],a2[i]] for i in range(len(xyz))]
    
    return dict(zip(serialn,volumes)), xyz1, xyz2

def calculateDM(coord1, coord2, maxCalcSize=5000, use_cuda=False):

    if use_cuda:
        darr = cuda.getDistanceMatrixSq3D(coord1, coord2, maxCalcSize)

    else:
        darr = npy.zeros([len(coord2), len(coord1)])
        for i in range(len(coord2)):
            darr[i][:]=(( coord1 - coord2[i]) ** 2).sum(axis=1)

    return darr
#
#def testASA(pdb, top):
#    import time
##    pdb = 'test/PEP_NAMD.pdb'
##    top = 'test/PEP_NAMD.top'
#
#    p = bi.PDBModel(pdb)
#    p.loadAmberTopology(top)
#    t0=time.time()
#    a, points = calculateASAamber(p, 1.4, n_sphere_point=150, getpoints=True, radii=True)
#    print "new:",time.time()-t0
##    t0=time.time()
##    a2, points2 = calculateASAamber_OLD(p, 1.4, n_sphere_point=150, getpoints=True, radii=True)
##    print "old:",time.time()-t0
#    printXYZ(points, 'testpoints.xyz')
##    printXYZ(points2, 'oldpoints.xyz')
#
#def testVolume():
#    pdb = 'test/PEP_NAMD.pdb'
#    top = 'test/PEP_NAMD.top'
#
#    p = bi.PDBModel(pdb)
#    p.loadAmberTopology(top)
#
#    volumes = getContactVolumes(p, 1.4, 4.9)
#
#    for serial,vol in volumes.iteritems():
#        print serial,': ',vol

#def buildGraph(points):
#
#    npoints = len(points)
#
#    # prepare the nodes
#    g = graph()
#    g.add_nodes(range(npoints))
#    [g.add_node_attribute(i, ('position',(points[i][0],points[i][1],points[i][2]))) for i in xrange(npoints)]
#
#    dm = calculateDM(points, points)
#
#    near = dm < 2.5**2
#    for i in range(npoints):
##        print near[i,:]
#        nh_ids = npy.where(near[i,:])[0]
##        print i, nh_ids
##        [g.add_edge(i,ni,wt=dm[i,ni]) for ni in nh_ids if ni != i]
#        for ni in nh_ids:
#            if ni != i:
##                print i,ni,dm[i,ni]
#                try:
#                    g.add_edge((i,ni),wt=npy.sqrt(dm[i,ni]))
#                except AdditionError:
#                    pass
#
#    h = euclidean.euclidean()
#    h.optimize(g)
#
#    return g, h
#    print g.node_neighbors
#    print dm[1,:]
#    dot = write(g)
#    gvv = gv.readstring(dot)
#    gv.layout(gvv, 'dot')
#    gv.render(gvv,'png','graph.png')
#    print h(30,33), dm[30,33]
#    tree, path = minmax.shortest_path(g, 0)
#    print tree, path
#    print minmax.heuristic_search(g, 0, 157,h), dm[0,157]
#
#def surfDist(visitlist,heur):
#    dist = [heur(visitlist[i],visitlist[i+1]) for i in range(len(visitlist)-1)]
##    print visitlist
##    dist = map(npy.sqrt,dist)
#    return reduce(lambda x,y: x+y, dist)

class surfaceGraph:
    def __init__(self, pdbModel, probe=1.4, n_sphere_points=30, edge_dist=2.5):
        self.pdb = pdbModel
        self.probe = probe
        self.sph_points = n_sphere_points
        self.edgedist = edge_dist
        self.visited = {}
        self._buildGraph()
        self._missingNodes()

    def _missingNodes(self):
        ok = False
        i = 0
        while not ok:
            tree, path = minmax.shortest_path(self.graph, i)
            # Check that node i has path to at least 95% of the nodes
            if len(path) > 0.95*len(self.allNodesCoord):
                self.missing = [j for j in self.graph.nodes() if j not in path.keys()]
                ok = True
            else:
                i += 1

    def nodeSurfdist(self, node1, node2, returnNodes=False):
        if node1 in self.visited.keys():
#            print self.visited
            d = self.visited[node1].get(node2)
            if d: return d
        elif node2 in self.visited.keys():
#            print self.visited
            d = self.visited[node2].get(node1)
            if d: return d
        else:
            tree, path = minmax.shortest_path(self.graph, node1)
            self.visited[node1] = path
            return path[node2]
        return None
    
    def _buildGraph(self):
        a, points = calculateASAamber(self.pdb,  self.probe, n_sphere_point=self.sph_points, getpoints=True)
        points = npy.array(points)
        npoints = len(points)
        
        # Remove those closer than 1 angstrom to each other
        dm = self._calculateDM(points, points)
        close = dm <= 1.5
        remove = []
        keep = []
        for i in xrange(npoints):
            if i not in remove:
                neighs = close[i,:]
                neighs[i] = False
                [remove.append(n) for n in npy.where(neighs)[0] if n not in keep]
                keep.append(i)
        points = npy.take(points, keep, axis=0)
        self.allNodesCoord = points
        npoints = len(points)
        
        # prepare the nodes
        g = graph()
        g.add_nodes(range(npoints))
        [g.add_node_attribute(i, ('position',(points[i][0],points[i][1],points[i][2]))) for i in xrange(npoints)]
        
        # add edges according to the distances between the points
        dm = self._calculateDM(points, points)
        self.nodesDM = dm

        # give high number to self distance instead of zero
        for i in range(len(dm)):
            dm[i,i] = npy.inf
        
        near = dm < self.edgedist
        for i in range(npoints):
            nh_ids = npy.unique(npy.where(near[i,:])[0])
            nn = len(nh_ids)
#            print i, "has",nn, "neighs"
            if nn > 0: 
                for ni in nh_ids:
                    if ni != i:
                        try:
                            g.add_edge((i,ni),wt=dm[i,ni])
                        except AdditionError:
#                            print "Could not add edge ",i,"to",ni
                            pass
            else:
#                print i,"has no neighs in edge distance"
                # No close nodes for self.edgedist chosen
                # will add one edge to the closer neighbour
                mydist = dm[i,:]
                neigh_i = npy.where(mydist == mydist.min())[0][0]
#                ok = False
                try:
                    g.add_edge((i,neigh_i),wt=dm[i,neigh_i])
#                    print "adding edge to ",neigh_i
                except:
#                    print "Could not add edge"
                    pass
#                print "adding edge to ",neigh_i
        
        # Do a final checking to prevent unconnected nodes
        # remove the unconected nodes at this point
#        rm_nodes = [node for node,neighs in g.node_neighbors.iteritems() if not neighs]
#        [g.del_node(ni) for ni in rm_nodes]
#        self.allNodesCoord = npy.delete(self.allNodesCoord, rm_nodes, 0)
        
        self.graph = g
    
    def _calculateDM(self, coord1, coord2):
        darr = npy.zeros([len(coord2), len(coord1)])
        if len(coord2) > 1:
            for i in range(len(coord2)):
                darr[i][:]=(( coord1 - coord2[i]) ** 2).sum(axis=1)
        else:
            darr = (( coord1 - coord2) ** 2).sum(axis=1)
            
        return npy.sqrt(darr)
        
    def getNodeXYZ(self, nodei):
        return self.allNodesCoord[nodei]
        
    def getCloserNode(self, coord, getDist=False):
        "Returns de closer node ID to the coordinate coord (numpy array)"
        ndist = npy.sqrt(((self.allNodesCoord - coord) ** 2).sum(axis=1))
        min = ndist.min()
        node = npy.where(ndist == min)[0][0]
        if getDist:
            return node, min
        else:
            return node

    def cartesianSurfDist(self, coord1, coord2, cut=2):
        "Returns distance between two points near the surface, following surface graph.\
        If the point is further than 'cut' distance from any node, it returns npy.inf."
        close1, d1 = self.getCloserNode(coord1, True)
        close2, d2 = self.getCloserNode(coord2, True)

        if d1 > cut or d2 > cut:
            return npy.inf
        if close1 == close2:
            return npy.sqrt(((coord2 - coord1)**2).sum())
        elif close1 in self.missing or close2 in self.missing:
            return npy.inf
#        try:
        nodeDist = self.nodeSurfdist(close1, close2)
#        d = nodeDist + d1 + d2
        d=nodeDist
#        except:
#            d = npy.nan
        return d
        
#def timedist(graph,ni,nj):
#    import time
#    t0=time.time()
#    print graph.nodeSurfdist(ni,nj)
#    print time.time() - t0
#
#def testGraph():
#    pdb = 'test/HEWL.pdb'
#    top = 'test/HEWL.top'
#    p = bi.PDBModel(pdb)
#    p.loadAmberTopology(top)
#    import time
#    t0=time.time()
#    pgraph = surfaceGraph(p)
#    print 'create graph:',time.time() - t0
#
#    return pgraph

if __name__ == "__main__":
#    testGraph()
#    testVolume()
#    import sys
#    testASA(sys.argv[1],sys.argv[2])
    pass