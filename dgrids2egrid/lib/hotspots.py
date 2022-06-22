#!/usr/bin/env python
# -*- coding: utf-8 -*-

from gridData import Grid
import numpy as np

def get_most_negative_points(g, threshold=-0.87):
    points = []
    for xi in range(len(g.grid)):
        for yi in range(len(g.grid)):
            for zi in range(len(g.grid)):
                if g.grid[xi][yi][zi] <= threshold:
                    coord = (g.origin[0]+g.delta[0]*xi, g.origin[1]+g.delta[1]*yi, g.origin[2]+g.delta[2]*zi)
                    points.append((coord,g.grid[xi][yi][zi]))
    return sorted(points, key=lambda x: x[1])

def distance(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(a-b)

def get_close_gridpoints(g, point_coord, reach = 2):
    xi = int(((point_coord[0] - g.origin[0])/g.delta[0])//1)
    yi = int(((point_coord[1] - g.origin[1])/g.delta[1])//1)
    zi = int(((point_coord[2] - g.origin[2])/g.delta[2])//1)

    points_around = [] # [((x,y,z), value), ...]
    deltas = range(-reach, reach+1)
    for dx in deltas:
        for dy in deltas:
            for dz in deltas:
                x_step = xi+dx
                y_step = yi+dy
                z_step = zi+dz
                if x_step >= g.grid.shape[0]: x_step = g.grid.shape[0]-1
                if y_step >= g.grid.shape[1]: y_step = g.grid.shape[1]-1
                if z_step >= g.grid.shape[2]: z_step = g.grid.shape[2]-1
                points_around.append((((g.origin[0]+g.delta[0]*(x_step)),
                                        (g.origin[1]+g.delta[1]*(y_step)),
                                        (g.origin[2]+g.delta[2]*(z_step))),
                                         g.grid[x_step][y_step][z_step]))
    return points_around

def filter_points_by_distance(ref, points, distance_cutoff = 2):
    new_list = [] # [((x,y,z), value, distance), ...]
    for point in points:
        d = distance(ref, point[0])
        if d <= distance_cutoff and point[1]<999:
            new_list.append((point[0], point[1], d))
            #print(point)
        elif d <= distance_cutoff and point[1]<= 999:
            pass
    return new_list

def isclose(hotspot, other_hotspots, distance_thr = 2):
    for other_hotspot in other_hotspots:
        d = distance(other_hotspot[0], hotspot[0])
        if d < 2:
            return True
    return False

def get_freenergy(g, hotspot, radius = 0.7):
    from math import e, log
    T = 300
    r = 0.0019872041
    coords = hotspot[0]
    points_sq = get_close_gridpoints(g, coords, reach=int((radius*2)+2))
    points = filter_points_by_distance(coords, points_sq, radius)
    #fenergy1 = -r*T*log(sum([e*(point[1]/-r*T) for point in points])/len(points))
    #fenergy = np.mean([x[1] for x in points])
    fenergy = get_boltzmannAVG([x[1] for x in points])
    #print(fenergy, fenergy2)
    return round(fenergy, 2)

def get_boltzmannAVG(data, T=300):
    from math import exp, log
    r = 0.0019872041

    #normalize
    reversedd = [-x for x in data]
    min_w = min(reversedd)
    norm_data = [x-min_w for x in reversedd]

    expWork = np.mean(list(map(lambda x: exp(x/(-r*T)), norm_data)))
    expAVG = T*r*log(expWork)-abs(min_w)
    return expAVG

def hotspots2pdb_single(file_handle, hotspots, probe_name = 'CT'):
    res_name = {'CT': 'HYD',
                 'OH': 'POL'}
    atom_name = {'CT': 'C',
                 'OH': 'O'} 
    atom_n = 0
    for hotspot in hotspots:
        atom_n += 1
        res_n = atom_n
        atom = atom_name[probe_name].center(4)
        res = res_name[probe_name].ljust(3)
        x = str('%8.3f'%hotspot[0][0]).rjust(8)
        y = str('%8.3f'%hotspot[0][1]).rjust(8)
        z = str('%8.3f'%hotspot[0][2]).rjust(8)
        occ = str('%6.2f'% float(1)).rjust(6)
        bfactor = str('%6.2f'%hotspot[1]).ljust(6)
        file_handle.write("%s%s %s %s %s%s    %s%s%s%s%s\n" %('ATOM'.ljust(6), str(atom_n).rjust(5),atom,res, 'A'.rjust(1), str(res_n).rjust(4), x, y, z, occ, bfactor ))

def hotspots2pdb(file_handle, hotspots_dict):
    res_name = {'CT': 'HYD',
                 'OH': 'POL'}
    atom_name = {'CT': 'C',
                 'OH': 'O'} 
    atom_n = 0
    for probe in hotspots_dict:
        hotspots = hotspots_dict[probe]
        for hotspot in hotspots:
            atom_n += 1
            res_n = atom_n
            atom = atom_name[probe].center(4)
            res = res_name[probe].ljust(3)
            x = str('%8.3f'%hotspot[0][0]).rjust(8)
            y = str('%8.3f'%hotspot[0][1]).rjust(8)
            z = str('%8.3f'%hotspot[0][2]).rjust(8)
            occ = str('%6.2f'% float(1)).rjust(6)
            bfactor = str('%6.2f'%hotspot[1]).ljust(6)
            file_handle.write("%s%s %s %s %s%s    %s%s%s%s%s\n" %('ATOM'.ljust(6), str(atom_n).rjust(5),atom,res, 'A'.rjust(1), str(res_n).rjust(4), x, y, z, occ, bfactor ))



def get_hotspots(g):
    good_hotspots = []
    hotspots = get_most_negative_points(g, threshold=-0.5)
    for hotspot in hotspots:
        if isclose(hotspot, good_hotspots): continue
        fenergy = get_freenergy(g, hotspot)
        good_hotspots.append((hotspot[0], fenergy))
    n_hotspots = int(((len(good_hotspots)*0.2)//1)+1)
    if n_hotspots > 400: n_hotspots = 400
    good_hotspots.sort(key=lambda x: x[1])
    return good_hotspots[:n_hotspots]

if __name__=='__main__':

    # params
    folder = 'PROBE_AVG'
    output_file = 'Hotspots_prueba.pdb'

    #load grids
    grids = {'CT' : Grid(folder+'/ETA_CT_DG0.dx'), 'OH': Grid(folder+'/ETA_OH_DG0.dx')}
    
    with open(output_file, 'w') as ofile:
        r = {}
        for probe in grids:
            g = grids[probe]
            r[probe] = get_hotspots(g)
        hotspots2pdb(ofile, r)

