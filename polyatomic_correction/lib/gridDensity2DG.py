import sys
import pyMDMix

if len(sys.argv) < 4:
	sys.exit("USAGE: gridDensity2DG.py inputgrid outputgrid expectednum [temp]\n grid input formats: dx, xplor, cns\ngrid output format: dx")

# Collect argument variables
ingrid = sys.argv[1]
outgrid = sys.argv[2]
expected = float(sys.argv[3])
try:
	temp = float(sys.argv[4])
except:
	temp = 300.

# Load grid and convert
grid = pyMDMix.Grid(ingrid)
dggrid = grid.count2DG(expected, temp)
grid.update(dggrid)
grid.writeDX(outgrid)

print "DONE"
