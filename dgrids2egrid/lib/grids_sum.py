#!/usr/bin/env python

from gridData import Grid
import sys

def gsum(A, B, C):
    D = A + B + C
    return D

if __name__=='__main__':

    grid_1 = sys.argv[1]
    grid_2 = sys.argv[2]
    grid_3 = sys.argv[3]
    out_grid = sys.argv[4]
    A = Grid(grid_1)
    B = Grid(grid_2)
    C = Grid(grid_3)
    D = gsum(A, B, C)
    D.export(out_grid)
