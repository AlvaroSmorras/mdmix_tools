#!/usr/bin/env python

from gridData import Grid
import numpy as np
from math import log 
def solvent_accessible_dots(g):
    non_accesible = 0
    for x in g.grid:
        for y in x:
            for z in y:
                if z == 0:
                    non_accesible += 1
    return (g.grid.shape[0]*g.grid.shape[1]*g.grid.shape[2])-non_accesible

def sum_counts(g):
    total_count = 0
    for x in g.grid:
        for y in x:
            for z in y:
                total_count += z
    return(total_count)

def calculate_egrid(dgrid, dgcorrection = 0):
    T = 300
    kb = 0.0019872041
    accesible_points = solvent_accessible_dots(dgrid)
    total_count = sum_counts(dgrid)
    # fenergy = -kb*T*log(point_count/(total_count/accesible_points))
    for x in range(len(dgrid.grid)):
        for y in range(len(dgrid.grid)):
            for z in range(len(dgrid.grid)):
                point_count = dgrid.grid[x][y][z]
                if point_count == 0:
                    fenergy = 999.999
                else:
                    fenergy = round(-kb*T*log(point_count/(total_count/accesible_points)), 3)
                dgrid.grid[x][y][z] = fenergy - dgcorrection
    return dgrid
if __name__=='__main__':
    folder = 'PROBE_AVG/'
    dpath = folder + 'ETA_WAT_dgrid.dx'
    epath = folder + 'ETA_WAT_DG.dx'
    dgrid = Grid(dpath)

    egrid = calculate_egrid(dgrid)
    egrid.export(epath)