#! /home/seorito/software/anaconda3/bin/python
#
#               surface_cluster.py
#
# This script calculates clusters of points using surface distances
# and graphs. It explores the points with lower energy and lets them
# grow a graph with the neighbouring points until a maximum of 2 levels of
# neighbours. Then calculates the energy of the graphs (suming up nodes contributions)

__author__="dalvarez"
__date__ ="$07-abr-2011 13:09:52$"

import sys
from pygraph.classes.graph import graph # python-graph-core
import Biskit as bi
import copy
import numpy as npy
from itertools import permutations
from os.path import exists

# globals
vdw_dict = {'C.3' : 1.88, 'C.2' : 1.76, 'C.1'   : 1.61, 'C.ar' : 1.88, 'C.cat' : 1.88,
            'N.3' : 1.64, 'N.2' : 1.64, 'N.1'   : 1.64, 'N.ar' : 1.64, 'N.am'  : 1.64, 'N.4' : 1.64, 'N.pl3' : 1.63,
            'O.3' : 1.46, 'O.2' : 1.42, 'O.co2' : 1.42, 'O.spc': 1.42, 'O.t3p' : 1.42,
            'S.3' : 1.782, 'S.2' : 1.77, 'S.O' : 1.77, 'S.O2' : 1.77, 'P' : 1.871,
            'F' : 1.56, 'Cl' : 1.735, 'Br' : 1.978, 'I' : 2.094, 'ZN' : 0.6, 'H': 0}

def get_sybyl_from_mol2(mol2):
    atom_flag = 0
    id_to_atomtype = {}
    with open(mol2) as mol2_fh:
        for line in mol2_fh:
            line = line.strip()
            if atom_flag:
                if line.startswith('@<TRIPOS>BOND'):
                    return id_to_atomtype
                line = line.split()
                id_to_atomtype[line[0]] = line[5]
                continue
            if line.startswith('@<TRIPOS>ATOM'):
                atom_flag = 1
            
def add_vdw(refPDB, mol2):
    atom_dict = get_sybyl_from_mol2(mol2)
    for atom in refPDB:
        atom['alternate'] = vdw_dict[atom_dict[str(atom['residue_number'])]]
def surfaceDM(refPDB, points, clean=False, graph=False, probe=1.2, edge_dist=3, n_sphere_points=50, verbose=False):
    """Returns a distance matrix for the list of points using surface distances.

    refPDB - is a Biskit.PDBModel of the reference structure over which surface
        distances will be calculated. AmberTopology must be loaded previously on the Model.

    points - is a npy.array (Nx3) with the coordinates to build the DM

    clean - if True, points not connected to the main graph (with inf distance) will be removed
        from the DM and a list of those indexes removed also returned

    graph - if True, the surfaceGraph instance is also returned as last argument

    SURFACE GRAPH Parameters (modify only in some special cases)

    probe - probe radius for ASA calculation. Used for surface graph building.
        Default = 1.2
    edge_dist - used in surface graph building. Distance to consider two vertex to be connected.
        Default = 3
    n_sphere_points - used in surface graph building. Number of points in the sphere for each atom.
        More points more dense, more accurate but much slower to calcualte. Surface graphs
        do not need such an accurate calculation.
        Default = 50

    It returns a npy.array (NxN)
    """

    from asa_calculations import surfaceGraph, printXYZ
    import time

    # Build surface graph
    t0 = time.time()
    if verbose: print( "Building surface graph...")
    sgraph = surfaceGraph(refPDB, probe=probe, edge_dist=edge_dist, n_sphere_points=n_sphere_points)
    if verbose: print( "Done in ",time.time() - t0)
    numpoints = len(points)
    surfDM = npy.zeros((numpoints,numpoints))

    printXYZ(sgraph.allNodesCoord, 'sgraph.xyz')

    # Calculate surface distances between all points
    t0=time.time()
    if verbose: print( "Calculating surface distance matrix...")
    N = len(points)
    for i, xyz1 in enumerate(points):
        surfDM[i,:] = [sgraph.cartesianSurfDist(xyz1, xyz2) for xyz2 in points]
        sys.stdout.write("%.2f"%(((i+1) / float(N)) * 100)+'% done\b\r');sys.stdout.flush()
    #print
    if verbose: print( "Done in ",time.time() - t0)

    # Return DM
    if clean:

        # Check how many inf values in the DM and what are the point with most of them
        # usually disconnected points from main graph should have small number of neighbours
        # and thus a number of inf values close to the total num of points
        if verbose: print( "Cleaning the DM removing disconected points")
        inf_mask = npy.isinf(surfDM)
        inf_sum = inf_mask.sum(axis=0)

        # The minimum number of inf will give us the number of disconnected points
        # All points with larger numinf will be removed. Here we assume that only one main graph
        # is generated and any smaller subgraph should be removed! TODO check if always correct.
        min_numinf = inf_sum.min()
        to_remove_ids = npy.where(inf_sum > min_numinf)[0]

        # Remove in both axis
        tmparr = npy.delete(surfDM, to_remove_ids, axis=0)
        surfDM = npy.delete(tmparr, to_remove_ids, axis=1)
        
        if verbose: print( "DONE")
        
        if graph: return surfDM, to_remove_ids, sgraph
        else: return surfDM, to_remove_ids
    
    else:
        if graph: return surfDM, sgraph
        else: return surfDM


def edges_graph(graph, dm, points_list=False, ranges_list=False, used=False):
    if points_list is False:
        points_list = graph.nodes()
    if ranges_list is False:
        ranges_list = range(len(graph.nodes()))

    if not used: used = []
    n = len(ranges_list)
    if n == 1:
        return graph
    else:
        rand = npy.random.randint(n)
        i = ranges_list[rand]
        used.append(i)
        d = dm[i,:]
        du = npy.delete(d, used)
        min = du.min()
        ranges_list = npy.delete(ranges_list, rand)

        closer_id = npy.where(d == min)[0][0]
        closer = points_list[closer_id]

        try:
            p = points_list[i]
            graph.add_edge((p,closer),wt=min)
        except:
            pass
        return edges_graph(graph, dm, points_list, ranges_list, used=used)

def getBestEdges(g, graphDM, tries=200):
    "Returns the graph with the best edges. Being that the one with the minimum sum of edges distances"
    import copy
    testg = []
    for i in range(tries):
        g_copy = copy.deepcopy(g)
        tmpg = edges_graph(g_copy, graphDM, g_copy.node_attr.keys())
        if isconnected(tmpg): testg.append(tmpg)
    allsums = [npy.sum([graph.edge_weight(i) for i in uniquePairList(graph.edges())]) for graph in testg]
    bestg = testg[npy.where(allsums == min(allsums))[0][0]]
    del(testg)
    return bestg

def saveGraph(g, name):
    dot = write(g)
    gvv = gv.readstring(dot)
    gv.layout(gvv,'dot')
    gv.render(gvv,'png',name)

def isconnected(g):
    from pygraph.algorithms.accessibility import connected_components as cc
    components = cc(g).values()
    uni_c = npy.unique(components)
    if len(uni_c) > 1: return False
    else: return True

def uniquePairList(pairlist):
    "Returns a list with the unique pairs in a list. Avoiding bidireccionalities."
    ue = []
    alle = []
    for ed in pairlist:
        if ed not in alle:
            ue.append(ed)
            alle.append(ed)
            alle.append((ed[1],ed[0]))
    return ue

def PDBstringFromGraph(g,chain='A',resid=1):
    pdbstr = 'ATOM  %5d%4s  %3s %1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f\n'
    type_elem = {'HYD':'C','DNA':'O','ACC':'F','DON':'N','POS':'S','NEG':'P', 'POL':'O'} # DNA is donor and acceptor?
    s = '' 
    for key,vals in g.node_attr.iteritems():
        vals = dict(vals)
        x,y,z = vals['position']
        energy = vals['energy']
        res = vals['type']
        at = type_elem[res]
        i = key
        s += pdbstr%(i,at,res,chain,resid,x,y,z,0.0,energy)
    return s

def PDBdescriptionFromGraph(g, chain='A'):
    pdbstr = 'HEAD -----------\nHEAD  GRAPH on chain %s:\nHEAD  DG: %.2f kcal/mol\nHEAD  Ki: %.4f nM\n'
    e = graphTotalEnergy(g)
    ki = DGtoKi(e, units='nM')
    if e < -15:
        pdbstr+='HEAD  WARNING, DG value too low, maybe clustering was greedy!\n'
    return pdbstr%(chain, e, ki)

def DGtoKi(dg, temp=300, units='M'):
    "Convert DG energy values to Ki values in the specified units\
    Default is M, other options: mM, uM, nM, pM"
    
    units_dict = {'M':1,'mM':1e3,'uM':1e6,'nM':1e9,'pM':1e12}
    R = 1.987/1000. # kcal/mol.K
    ki = npy.exp(dg/(R*temp))
    
    return ki*units_dict[units]

def writeAllClusters(g_list, pdbname):
    hs = ''
    ms = ''
    import string
    ascii = string.ascii_uppercase

    # Sort g_list by lower to higher energies
    en = [graphTotalEnergy(g) for g in g_list]
    to_sort = zip(en, g_list)
    to_sort.sort()
    
    for i,p in enumerate(to_sort):
        if i > 25: k = i - 26
        else: k = i
        chain = ascii[k]
        hs += PDBdescriptionFromGraph(p[1],chain=chain)
        ms += PDBstringFromGraph(p[1], chain=chain, resid=i+1)
    f=open(pdbname,'w')
    f.write(hs+ms)
    f.close()

def graphTotalEnergy(g):
    e = npy.array([dict(vals)['energy'] for vals in g.node_attr.values()])
    return e.sum()

def euclideanDist(p1, p2):
    return npy.sqrt(((p2 - p1)**2).sum())

def graphXYZ(g):
    gxyz = npy.zeros((len(g),3))
    for i,node in enumerate(g.nodes()):
        gxyz[i,:] = dict(g.node_attributes(node))['position']
    return gxyz

def eDistEdge(g,edge):
    n1,n2=edge
    "Return the euclidean distance between a pair of nodes"
    pos1 = npy.array(dict(g.node_attributes(n1))['position'])
    pos2 = npy.array(dict(g.node_attributes(n2))['position'])
    return euclideanDist(pos1, pos2)

def graphVolume(g):
    """
    Calculates a very approximate volume of a molecule touching all the nodes.
    Each edge will be considered as an atom chain. One atom every 1.5 Anstroms.
    Radii of 1.9 A per atom. V = (4/3)*pi*r^3 = 28.73 A^3.
    """
    ATOM_VOL = 28.73
    g_atoms = len(g.nodes())
    g_edges = uniquePairList(g.edges())
#    edge_dist = npy.array([g.edge_weight(e) for e in g_edges])
    edge_dist = npy.array([eDistEdge(g,e) for e in g_edges])
    linking_atoms = npy.floor(edge_dist / 1.5).astype(int) - 1
    linking_atoms = linking_atoms.sum()
    volume = (linking_atoms + g_atoms)*ATOM_VOL
    return volume

def distBetweenGraphs(g1,g2,dm, nodes=False):
    "Returns the shortest distance between two graphs nodes using distance matrix 'dm' \
    containing all initial coordinates before clustering (surfDM).If nodes is TRUE, return\
    also closer node ids. "
    if g1 == g2:
        return 0
    n1 = g1.nodes()
    n2 = g2.nodes()
    A = npy.take(dm, n1, axis=0)
    A = npy.take(A, n2, axis=1)
    d = A.min()
    if nodes:
        x,y = npy.where(A == d)
        pair = (n1[x[0]], n2[y[0]])
        return d, pair
    else:
        return d

def connect2Graphs(g1,g2,n1,n2,wt=False):
    "Returns a new graph object combining all graph1 and graph2 by connecting them\
    with an edge from n1 node in g1 to n2 node in g2. give a wt to the new edge aswell, or it will be zero"
    from pygraph.classes.graph import graph
    g = graph()
    # crete new graph with existing nodes and edges
    for gr in (g1,g2):
        for node, attr in gr.node_attr.iteritems():
            g.add_node(node, attrs=attr)
        for edge in uniquePairList(gr.edges()):
            w = gr.edge_weight(edge)
            g.add_edge(edge, wt=w)
    # Add new edge
    if not wt: wt=1
    g.add_edge((n1,n2),wt=wt)
    return g

def getXyzTypesEnergies(pointsPDB_file):
    
    pointsPDB = bi.PDBModel(pointsPDB_file)
    m = npy.bitwise_not(pointsPDB.maskFrom('residue_name','WAT'))
    points = pointsPDB.compress(m)
    xyz = points.xyz
    types = points['residue_name']
    energies = points['temperature_factor']

    return xyz, types, energies

def buildsurfDM(refPDB, xyz, surfDMfile, removedfile, FORCE=False):
    # Calculate surface distance matrix
    # Check if pickles already exists to reduce calcualtions
    # FORCE can force re-calculation independently of the existence of the files
    if exists(surfDMfile) and exists(removedfile) and not FORCE:
        # Load Precalculated matrix
        print( "Loading precalculated surface distance matrix...")
        surfDM = bi.tools.load(surfDMfile)
        removed_ids = bi.tools.load(removedfile)
    else:
        # Calculate matrix and store as pickle
        surfDM, removed_ids = surfaceDM(refPDB, xyz, clean=True, verbose=False)
        bi.tools.dump(surfDM,surfDMfile)
        bi.tools.dump(removed_ids, removedfile)

    return surfDM, removed_ids

def clusterPoints(pdbfile, topfile, pointsfile, name, neigh_cut=6, adjacent_cut=8, POINTSPERGRAPH=2, FORCE=False):

    # Set file names
    surfDMfile = name+"_surfDM.pick"
    removedfile = name+"_removed.pick"
    
    # Load files and prepare instances
    print( "Loading files")
    refPDB = bi.PDBModel(pdbfile)
    #refPDB.loadAmberTopology(topfile)
    add_vdw(refPDB, topfile)

    # Get all information about the points to cluster
    xyz, types, energies = getXyzTypesEnergies(pointsPDB_file)

    # Calculate surfce distance matrix and obtain list of discarded points
    # (those far from any surface node) or innaccessible from the surface (inner closed cavities).
    surfDM, removed_ids = buildsurfDM(refPDB, xyz, surfDMfile, removedfile, FORCE=FORCE)

    # Remove ids also from types, energies, xyz and ids
    DM = surfDM[:,:]
    types = npy.delete(types, removed_ids)
    energies = npy.delete(energies, removed_ids)
    xyz = npy.delete(xyz, removed_ids, axis=0)
    ids = npy.arange(len(xyz))

    # Loop over all points and get maximum number of sub-clusters
    pocket_graphs = []
    while len(xyz):
        min_id = npy.where(energies == energies.min())[0][0]
        min_neigh_ids = npy.where(DM[min_id,] <= neigh_cut)[0]
        # If no neighbours, remove from all arrays
        # length will be 1 (itself)
        if len(min_neigh_ids) == 1:
            types = npy.delete(types, min_id)
            energies = npy.delete(energies, min_id)
            xyz = npy.delete(xyz, min_id, axis=0)
            ids = npy.delete(ids, min_id)
            DM = npy.delete(DM, min_id, axis=0)
            DM = npy.delete(DM, min_id, axis=1)
            continue
        # Check neighbours of neighbours
        subDM = npy.take(DM, min_neigh_ids, axis=0)
        subzone_ids = npy.where(npy.sum(subDM <= neigh_cut, axis=0) != 0)[0]
        #permisive_zone_ids = npy.where(npy.sum(subDM <= permisive_cut, axis=0) != 0)[0]
        # Take all info for the selected points
        # build a sub distance matrix to later build a graph
        graph_egs = energies[subzone_ids]
        graph_types = types[subzone_ids]
        graph_xyz = xyz[subzone_ids]
        graph_ids = ids[subzone_ids]
        graphDM = npy.take(DM, subzone_ids, axis=0)
        graphDM = npy.take(graphDM, subzone_ids, axis=1)
        graphDM[graphDM==0] = npy.inf   # put diagonal to infinity
        graphDM = graphDM.astype("float32")
        # remove from pool of hotspots the selected ones
        types = npy.delete(types, subzone_ids)
        energies = npy.delete(energies, subzone_ids)
        xyz = npy.delete(xyz, subzone_ids, axis=0)
        ids = npy.delete(ids, subzone_ids)
        DM = npy.delete(DM, subzone_ids, axis=0)
        DM = npy.delete(DM, subzone_ids, axis=1)
        # create graphs with more than 3 nodes
        n_elem = len(graph_egs)
        if n_elem >= POINTSPERGRAPH:
            # create graph and add all nodes and attributes
            # calculate the best edges between nodes
            g = graph()
            g.add_nodes(graph_ids)
            [g.add_node_attribute(n, ('position',(graph_xyz[i,0],graph_xyz[i,1],graph_xyz[i,2]))) for i,n in enumerate(graph_ids)]
            [g.add_node_attribute(n, ('energy',graph_egs[i])) for i,n in enumerate(graph_ids)]
            [g.add_node_attribute(n, ('type',graph_types[i])) for i,n in enumerate(graph_ids)]
            g = getBestEdges(g, graphDM)
            # append graph to the graph list
            pocket_graphs.append(g)

    # Try to join adjacent clusters
    # Find clusters closer than adjacent_cut distance
    ngraphs = len(pocket_graphs)
    print( "Found "+ str(ngraphs) + " subgraphs")
    print( "Joining those at "+str(adjacent_cut)+'A')
    pocket_dist = npy.zeros((ngraphs,ngraphs))
    for i,p in enumerate(pocket_graphs):
        pocket_dist[i,:] = [distBetweenGraphs(p,pi,surfDM) for pi in pocket_graphs]

    pocket_dist[pocket_dist == 0] = npy.inf
    c1, c2 = npy.where(pocket_dist <= adjacent_cut)
    graphpairs = uniquePairList(zip(c1,c2))
    print('\tjoining: ',graphpairs)

    # Find closer nodes between those graphs
    nodepairs = []
    for pair in graphpairs:
        x,y=pair
        d, n = distBetweenGraphs(pocket_graphs[x],pocket_graphs[y],surfDM,nodes=True)
        nodepairs.append(n)

    pairs_nodes = dict(zip(graphpairs,nodepairs))

    # Detect if same graph joins more than 1 different graph
    uni_graph_list = []
    [[uni_graph_list.append(i) for i in pair if i not in uni_graph_list] for pair in graphpairs]
    if len(uni_graph_list) < len(graphpairs)*2:
        # some of the graphs is repeated
        # find which
        rep = npy.zeros((len(uni_graph_list),len(graphpairs))).astype(bool)
        for i,u in enumerate(uni_graph_list):
            rep[i,:] = [u in pair for pair in graphpairs]
        repeated_ids = npy.where(rep.sum(axis=1) > 1)[0]
        repeated_graphs = npy.take(uni_graph_list, repeated_ids)
    else:
        repeated_graphs = []
    if len(repeated_graphs) > 0:
        print( "Found multiple connected graphs: ",repeated_graphs)

    # Check that a there is no more than 3 connections per graph
    if len(repeated_graphs) > 1:
        for perm in permutations(repeated_graphs, 2):
            if perm in graphpairs:
                print( "WARNING! More than 3 graphs connected in a row! Pair:",perm)

    # Connect pair of graphs which are not in repeated_graphs
    # make a list with the ones repeated
    connect_multiple = []
    new_connected = []
    for pair in graphpairs:
        if pair[0] in repeated_graphs or pair[1] in repeated_graphs:
            connect_multiple.append(pair)
        else:
            print( "Connecting pair:",pair)
            new_connected.append(connect2Graphs(pocket_graphs[pair[0]],pocket_graphs[pair[1]],pairs_nodes[pair][0],pairs_nodes[pair][1]))

    # Now connect the ones with multiple neighbours
    if connect_multiple:
        for m in repeated_graphs:
            myg = copy.deepcopy(pocket_graphs[m])
            for p in connect_multiple:
                if m in p:
                    print( "Multiple connecting: ",p)
                    p = npy.array(p)
                    other = p[p != m][0]
                    # get nodes in right order
                    if (m,other) in pairs_nodes.keys():
                        nodes = pairs_nodes[(m,other)]
                    else:
                        nodes = pairs_nodes[(other,m)]
                        nodes = (nodes[1], nodes[0])
                    myg = connect2Graphs(myg, pocket_graphs[other],nodes[0],nodes[1])
            new_connected.append(myg)

    # append the new generated clusters and remove the previous ones
    pocket_graphs += new_connected
    unique_graph_ids = []
    for pair in graphpairs:
        for i in pair:
            if i not in unique_graph_ids: unique_graph_ids.append(i)
    unique_graph_ids.sort(reverse=True)

    [pocket_graphs.pop(i) for i in unique_graph_ids]

    rm_repeated_nodes(pocket_graphs)
    return pocket_graphs
def rm_repeated_nodes(pocket_graphs):
    # repeated nodes in graphs
    nodes_in_pocket = {}
    all_nodes = set()
    for i, pocket in enumerate(pocket_graphs):
    	nodes_in_pocket[i] = set()
    	for a in pocket.node_attr.iteritems():
    		nodes_in_pocket[i].add(a[0])
    		all_nodes.add(a[0])

    
    all_nodes = list(all_nodes)
    node_pocket = []
    for node in all_nodes:
    	pocketss = []
    	for pocket in nodes_in_pocket:
    		if node in nodes_in_pocket[pocket]:
    			pocketss.append(pocket)
    	node_pocket.append((node, pocketss))
    
    node_pocket.sort(key= lambda x:x[1])
    #for x in node_pocket: print(x)

    
    for pocket in nodes_in_pocket:
    	for other_pocket in nodes_in_pocket:
    		if pocket == other_pocket: break
    		for node in nodes_in_pocket[pocket]:
    			for n in pocket_graphs[other_pocket]:
    				if n == node:
    					pocket_graphs[other_pocket].del_node(node)

def defineOptParser():
    parser = OptionParser()
    
    # MANDATORY
    parser.add_option("-n","--name", dest='name',
            help="Project name. Used as prefix for output file.")
    parser.add_option("-r","--ref", dest="ref",
            help="Reference structure PDB file name.")
    parser.add_option("-t","--top",dest="top",
            help="Mol2 file of the reference PDB, from which the atomtypes are extracted.")
    parser.add_option("-p","--points",dest="points",
            help="Points to cluster")

    # OPTIONAL
    group = OptionGroup(parser, "Optional args",
                    "Use them if you know what you are doing :)")

    group.add_option("-f","--force",dest="force",action="store_true",default=False,
            help="Force recalculation of surface distance matrix. Default False.")
    group.add_option("-i","--neighc",dest="neigh",type="float",default=6,
            help="Neighbour cutoff distance for node connection in subgraphs. Surface distance. Default=6 A")
    group.add_option("-j","--joinc",dest="join",type="float",default=8,
            help="Neighbour cutoff distance for subgraph joining. Surface distance. Default=8 A")
    group.add_option("-m","--min",dest="min",type="int",default=3,
            help="Minimum number of nodes to consider a subgraph. Default = 3 nodes.")

    parser.add_option_group(group)
    
    return parser

def parseArgs(parser, args):

    (options, uargs) = parser.parse_args(args)

    return options, uargs

def test():
    pdbfile = 'test/hsp90_reference.pdb'
    topfile = 'test/hsp90.prmtop'
    pointsPDB_file = 'test/hsp90_drug_allnoions_hotspots.pdb'
    name = 'HSP90'
    neigh_cut = 6
    adjacent_cut = 8
    POINTSPERGRAPH = 3
    FORCE = False
    pocket_graphs = clusterPoints(pdbfile, topfile, pointsPDB_file, name,
        neigh_cut=neigh_cut, adjacent_cut=adjacent_cut, POINTSPERGRAPH=POINTSPERGRAPH, FORCE=FORCE)
if __name__ == "__main__":
    #test()
    if 1:
        from optparse import OptionParser, OptionGroup

        parser = defineOptParser()

        if len(sys.argv) < 4:
            print( '-'*40)
            parser.print_help()
            sys.exit(1)

        opts, args = parseArgs(parser, sys.argv[1:])
        
        # User arguments
        name = opts.name
        pdbfile = opts.ref
        topfile = opts.top
        pointsPDB_file = opts.points
        
        # Parameters
        neigh_cut = opts.neigh
        adjacent_cut = opts.join
        POINTSPERGRAPH = opts.min
        FORCE = opts.force
        
        # RUN CLUSTERING
        pocket_graphs = clusterPoints(pdbfile, topfile, pointsPDB_file, name,
            neigh_cut=neigh_cut, adjacent_cut=adjacent_cut, POINTSPERGRAPH=POINTSPERGRAPH, FORCE=FORCE)
        
        # SAVE RESULTS
        outname = name+'_clusts.pdb'
        print( "Saving results: ",outname)
        writeAllClusters(pocket_graphs, outname)
