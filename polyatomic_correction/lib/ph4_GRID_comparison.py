#! /usr/bin/env python
"""
    Compare MPI GRID probes Grids with pharmacophore points saved as ph4 files.


	USAGE: ph4_GRID_compare.py inputph4.ph4 ph4ranking.txt path/to/grids.xplor outputprefix [onlycenter]

    Will use two metrics:
	1) Distance based metric from hotspot center coordinate
	2) Energy metric

	1) DISTANCE BASED METRIC
		Calculate hotspots from the input grids for each probe type and compare distances with
		Pharmacophore centers and radius. 
		
		If grid probe type matches ph4 type:
			+ If nearest hotspot inside ph4: +1
			+ If outside: -0.5
		If grid probe type differes from ph4 type:
			+ If nearest hotspot inside: -1
			+ If nearest hp outside: +0.5
	
	2) ENERGY METRIC
		Calculate a mean energy of the values enclosed by the ph4 (center + radius sphere).
		If types do not match: Invert energy value (as negative energy in non matching ph4 is negative
		and positive energies in non matching ph4 is positive). If types match leave energy as it is (negative energy is good
		while positive energy is bad).

    Will output 2 files:
	- outprefix_dist.csv 		with distance metric results
	- outprefix_ene.csv		with energy metric results
"""
import sys
import numpy as np
import Biskit as bi
import glob
import scipy.spatial as spat
import os
import matplotlib.pyplot as plt

import MDMix.GridsManager as G
import MDMix.HotSpotsManager as HM

###DICTIONARIES
solvent_types = {''}
#solvent_types = {'1PR':{'CT':['Hyd', 'Aro', 'PiR'],'O':['Acc','Don', 'Acc&Don','Don&Acc']},
#'ANT':{'C':['Hyd','Aro','PiR'],'N':['Acc']},
#'CLE':{'CL':['Hyd','Aro','PiR'],'CT':['Hyd','Aro','PiR'],'OH':['Acc','Don', 'Acc&Don','Don&Acc']},
#'ETA':{'CT':['Hyd','Aro','PiR'],'OH':['Acc','Don', 'Acc&Don','Don&Acc']},
#'IMZ':{'CA':['Aro','Hyd','PiR'],'N':['Acc'],'NH':['Don']},
#'ISX':{'CA':['Aro','Hyd','PiR'],'N':['Acc'],'O':['Acc']},
#'MAM':{'CT':['Hyd','Aro','PiR'],'N':['Don'],'O':['Acc']},
#'MSU':{'CA':['Hyd'],'N':['Don'],'O':['Acc']},
#'PYR':{'CA':['Aro','Hyd','PiR'],'N':['Acc']},
#'PYZ':{'CA':['Aro','Hyd','PiR'],'N':['Acc'],'NH':['Don']},
#'TFE':{'CF':['Hyd','Aro','PiR'],'O':['Acc','Don','Acc&Don','Don&Acc']},
#'ION':{'POS':['Ani'], 'NEG':['Cat']},
#}
#
### RegExp build to detect probes in filenames
import re
regexp =[]
for solv, probes in solvent_types.iteritems():
	regexp.extend([re.compile(".(%s_%s)."%(solv, p)) for p in probes.keys()])

def fetchprobe(s):
	"Return probe matched in string s"
	for reg in regexp:
		match = reg.search(s)
		if match: return match.groups()[0]
	return False


if len(sys.argv) < 5:
	sys.exit("USAGE: ph4_mdmix_superposition.py inputph4.ph4 ph4ranking.txt path/to/grids.dx outputprefix [onlycenter]")

### READ PH4 FILE 
inputph4 = sys.argv[1] # Name of the pharmacophore input (e.g. 'all_1st.ph4')
ph4rank = sys.argv[2]
inputpath = sys.argv[3] # Path of all grid files (e.g. 'HS_avg_pdb/*.dx')
outprefix = sys.argv[4]  # e.g. 'Final_results'
try:
	o = sys.argv[5]
	if o.lower() == 'onlycenter': onlycenter=True
	else: onlycenter=False
except:
	onlycenter=False


def readph4(inputph4):
	"Read ph4 formated file and return a list with all pharmacophoric point information"
	f = open(inputph4,'r') # Read the inputph4 file
	# Skip lines until #feature is found
	phs = []
	while True:
		l = f.readline()
		if '#feature' in l: break

	# Read all lines until # is found
	phlines = []
	while True:
		l = f.readline()
		if '#' in l: break
		phlines.append(l)

	f.close()
	phline = ' '.join(phlines) # join all lines with blank space to construct a single sentence
	phinfo = phline.split()

	# Function to split list l in chunks of n elements
	def chunks(l, n):
	    """ Yield successive n-sized chunks from l.
	    """
	    return [l[i:i+n] for i in xrange(0, len(l), n)]

	# Each pharmacophore consists of 9 elements of which 
	# 1st is Type, 3,4,5 are xyz and 6 is radius
	return chunks(phinfo, 9) #split list of elements in groups of 9
	

def readph4rank(f):
	out = []
	lines = open(f, 'r').readlines()
	[out.extend(l.strip().split()) for l in lines]
	return map(float, out)

#Read PHARMACOPHORE and save solvent types, xyz and radius in lists
phs = readph4(inputph4)
phs_rank = readph4rank(ph4rank)

# SORT ph4s according to rank
order = np.argsort(phs_rank)
phs = np.array(phs)[order].tolist()[::-1]
phs_rank = np.array(phs_rank)
phs_rank.sort()
phs_rank = phs_rank[::-1]

types = []
coords = []
ph4radius = []
for pharmacophore in phs: 
	phtype = pharmacophore[0].replace('$m','&')
	coord = map(float, pharmacophore[2:5]) #convert list of strings in list of numbers
	rad = float(pharmacophore[5]) #convert string to number
        types.append(phtype)   
        coords.append(coord)
	ph4radius.append(rad)

xyz_ph4 = np.array(coords) #convert list in numpy array
types_ph4 = [k.split('|') for k in types] #split by | 
# Obtain a unique list of types
allph4_types = []
[allph4_types.extend(k.split('|')) for k in types]
allph4_types = np.unique(allph4_types)
print 'Pharmacophore DONE'

########### Load grids, convert to Hotspots and evaluate energy
path = glob.glob(inputpath+'/*.dx')  #Read allpdbs in inputpath (e.g. HS_avg_pdb/*.pdb)
if not path:
	sys.exit("No files with extension *.dx found in folder %s"%inputpath)
results_distance = {}     #Dict of PROBE & Score from every Pharmacophore vs HSpoint comparison (e.g. results={'Probe':[score]})
results_energy = {}	 # Save energy inside ph4 point
tolerance=0.5    #Radius tolerance created to catch more points. 
RT = 0.0019859*300 # kcal/mol

for gridfile in path:
	# Identify probe from filename and get assigned types
	probe = fetchprobe(os.path.basename(gridfile))
	if not probe:
		print "File probe not identified. Skipping: %s"%gridfile 
		continue # skip if probe not identified

	Solv, AT = probe.split('_')  #Solvent types
        probe_types = solvent_types[Solv][AT] #properties (e.g. ['Acc'], ['Don'], ['Acc','Don'], ['Hyd'], ['Aro','Hyd'])
	if not	set(allph4_types) & set(probe_types):
		print "Probe %s does not have representation in pharmacophore file. Skipping..."%probe
		continue

	print "Analyzing grid %s, probe %s"%(gridfile, probe)
	# Load the grid file 
	g = G.Grid(gridfile)

	# Identify hotspots near each ph4 point by distance
	hmanage = HM.HotSpotsManager()
        hmanage.createHotSpotSetFromGrid(g, method='cutoff', name='')
	xyz_hs = []
	if onlycenter: [xyz_hs.extend([hotspot.coord]) for hotspot in hmanage.hotspotSet.hotspots]
	else: [xyz_hs.extend(hotspot.coordList) for hotspot in hmanage.hotspotSet.hotspots]
	xyz_hs = np.array(xyz_hs)
	md = spat.distance_matrix(xyz_ph4,xyz_hs) #A distance matrix measuring the distance between pharmacophore and HS coord. 
	partial_results=[]
	for q, score in enumerate(md): #take the minimum distance from every array.
		mdmin=score.min()
		ph4type=types_ph4[q]
		R=ph4radius[q]
		if set(ph4type).intersection(set(probe_types)):       # fc4 == ppt_hs
	       		if mdmin < (R + tolerance):
               			puntuation = '1'
			else:
               			puntuation = '-0.5'
		else:
			if mdmin < (R + tolerance):
				puntuation = '-1'
			else:       
				puntuation = '0.5'
		partial_results.append(puntuation)
	results_distance[probe]=partial_results	

	# Get grid energy value inside each ph4 point
	meanenergy_list = []
	for i, coord in enumerate(xyz_ph4):
		ph4type = types_ph4[i]
		if set(ph4type) & set(probe_types): match = True
		else: match=False
		rad = ph4radius[i]
		evals = g.getSphereValues(coord, rad)
		e=-RT*np.log(np.exp(evals / -RT).mean())
		# If types match, leave energy as it is
		# negative values are good and positive values are bar
		# if types do not match, swap signs to make negative positive (bad behaviour, good energy in non corresponding ph4)
		# or bad energy in non corresponding ph4 is a good indicator.
		if not match: e*=-1 # If incorrect type and energy is positive (good behaviour, invert energy value)
		meanenergy_list.append(e)
		
	results_energy[probe] = meanenergy_list
	

## PROCESS RESULTS
## Create a Table in .csv format separating strings by ';' and lines by \n (means new line)
outdistance = []
types = ['%s_%.2f'%(t,phs_rank[i]) for i,t in enumerate(types)]
header = ['Probes|Ph4']+types
outdistance.append(';'.join(header))
for key, values in results_distance.iteritems():
	line = [key]+values
	outdistance.append(';'.join(line))
outdistance = '\n'.join(outdistance)

# Build same with energy values
outenergy = [';'.join(header)]
for key, values in results_energy.iteritems():
	vals = ['%.3f'%v for v in values]
	line = [key]+vals
	outenergy.append(';'.join(line))
outenergy = '\n'.join(outenergy)

# SAVE RESULT TABLES
outname = outprefix+'_dist.csv'
print "Writting distance results CSV table:",outname
out_results = open(outname, 'w')
out_results.write(outdistance)
out_results.close()

outname = outprefix+'_ene.csv'
print "Writting energy results CSV table:",outname
out_results = open(outname, 'w')
out_results.write(outenergy)
out_results.close()
print "DONE"
#####################


################ PLOTTING PART
#### Heatmap plot of each of the results
#
#
## Normalize colorbar to ensure zero is white
#from matplotlib.colors import Normalize
#
#class MidpointNormalize(Normalize):
#    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#        self.midpoint = midpoint
#        Normalize.__init__(self, vmin, vmax, clip)
#
#    def __call__(self, value, clip=None):
#        # I'm ignoring masked values and all kinds of edge cases to make a
#        # simple example...
#        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#        return np.ma.masked_array(np.interp(value, x, y))
#norm = MidpointNormalize(midpoint=0)
#
#
#########################
## Load data into numpy arrays and append total score column
#dist_data = np.array(results_distance.values(),dtype='float')
#ene_data = np.array(results_energy.values(),dtype='float')
## Normalize energy data
##ene_data = (ene_data - ene_data.mean(axis=1).reshape(len(ene_data),1))/(ene_data.max(axis=1)-ene_data.min(axis=1)).reshape(len(ene_data),1)
#ene_totals = ene_data.sum(axis=1).reshape((len(ene_data),1))
#dist_totals = dist_data.sum(axis=1).reshape((len(dist_data),1))
#ene_data = np.hstack((ene_data, ene_totals))
#dist_data = np.hstack((dist_data, dist_totals))
#
#order = np.argsort(ene_data[:,-1])
#ene_data=ene_data[order]
#dist_data=dist_data[order]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(111) 
## Turn off axis lines and ticks of the big subplot
#ax.spines['top'].set_color('none')
#ax.spines['bottom'].set_color('none')
#ax.spines['left'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.get_xaxis().set_visible(False)
#ax.get_yaxis().set_visible(False)
#ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#
#def addheatmap(mainfig, numsubfigs, subfig, data, labelsx, labelsy, title):
#	ax1 = mainfig.add_subplot(numsubfigs,1,subfig)
#	# Get minimum and maximum values without total score
#	minval = data[:,:-1].min()
#	maxval = data[:,:-1].max()
#	if subfig==1: 
#		cm = 'PuOr'
#		p = ax1.pcolor(data, cmap=cm, alpha=0.8, edgecolors='w', linewidth=0.1, vmin=minval, vmax=maxval)	
#	else: 
#		cm = 'PuOr_r'
#		p = ax1.pcolor(data, cmap=cm, norm=norm, alpha=0.8, edgecolors='w', linewidth=0.1, vmin=minval, vmax=maxval)	
#	# put the major ticks at the middle of each cell
#	ax1.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
#	ax1.set_yticks(np.arange(data.shape[0])+0.5, minor=False)
#
#	# want a more natural, table-like display
#	ax1.invert_yaxis()
#	ax1.xaxis.tick_top()
#	ax1.set_frame_on(False)
#
#	#ax1.set_title(title)
#	if subfig == 1: ax1.set_xticklabels(labelsx, minor=False, rotation=90, fontsize=8)
#	else: ax1.set_xticklabels([], minor=False)
#	ax1.set_yticklabels(labelsy, minor=False, fontsize=8)
#	
#	# Remove ticks everywhere
#	for t in ax1.xaxis.get_major_ticks():
#	    t.tick1On = False
#	    t.tick2On = False
#	for t in ax1.yaxis.get_major_ticks():
#	    t.tick1On = False
#	    t.tick2On = False
#
#	# Add value to Total Score column
#	for y in range(data.shape[0]):
#    		x = data.shape[1]-1
#       		plt.text(x + 0.5, y + 0.5, '%.1f' % data[y, x],
#               	horizontalalignment='center',
#               	verticalalignment='center', fontsize=8)
#
#	cbar=plt.colorbar(p)	
#	cbar.ax.tick_params(labelsize=8)
#
## Add total score label to types
#types += ['TotalScore']
#distlabels = np.array(results_distance.keys())[order].tolist()
#enelabels = np.array(results_energy.keys())[order].tolist()
#addheatmap(fig, 2, 1, dist_data, types, distlabels, 'Distance scoring')
#addheatmap(fig, 2, 2, ene_data, types, enelabels, 'Energy instide ph4')
##ax.set_title('Scoring MDMix probes by pharmacophore overlapping')
#plt.tight_layout()
#plt.savefig(outprefix+'_heatmaps.pdf')
##plt.show()
