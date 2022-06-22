#! /usr/bin/env python
#
#  PLOT RESULTS FROM ph4_mdmix_superposition.py
#
#	Each probe is given a Score: each value is multiplied by the relevance of the ph4 point then all partial score sumed up
#
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

if len(sys.argv) < 4:
	sys.exit("USAGE: plot_ph4_mdmix_comparison.py ph4results_distance.csv ph4results_energy.csv outputplot.pdf") 

distancecsvf = sys.argv[1]
energycsvf = sys.argv[2]
outname = sys.argv[3]

# Transform csv to dictionary
results_distance = {}
results_energy = {}
dl = open(distancecsvf,'r').readlines()
phtypes = dl[0].split(';')[1:]
for l in dl[1:]:
	s = l.split(';')
	results_distance[s[0]] = np.array(s[1:],dtype=float)

el = open(energycsvf,'r').readlines()
for l in el[1:]:
        s = l.split(';')
        results_energy[s[0]] = np.array(s[1:],dtype=float)

############### PLOTTING PART
### Heatmap plot of each of the results
# Normalize colorbar to ensure zero is white
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
norm = MidpointNormalize(midpoint=0)


########################
# Load data into numpy arrays and append total score column
dist_data = np.array(results_distance.values(),dtype='float')
ene_data = np.array(results_energy.values(),dtype='float')

# Normalize energy data
#ene_data = (ene_data - ene_data.mean(axis=1).reshape(len(ene_data),1))/(ene_data.max(axis=1)-ene_data.min(axis=1)).reshape(len(ene_data),1)

# CALCULATE SCORES
# Extract relevance from ph4 titles
relevance = np.array([p.split('_')[1] for p in phtypes], dtype=np.float)
ene_score = (ene_data*relevance).sum(axis=1).reshape((len(ene_data),1))
dist_score = (dist_data*relevance).sum(axis=1).reshape((len(dist_data),1))
ene_data = np.hstack((ene_data, ene_score))
dist_data = np.hstack((dist_data, dist_score))

eorder = np.argsort(ene_data[:,-1])
ene_data=ene_data[eorder]
dorder = np.argsort(dist_data[:,-1])[::-1]
dist_data=dist_data[dorder]

fig = plt.figure()
ax = fig.add_subplot(111)
# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

def addheatmap(mainfig, numsubfigs, subfig, data, labelsx, labelsy, title):
        ax1 = mainfig.add_subplot(numsubfigs,1,subfig)
        # Get minimum and maximum values without total score
        minval = data[:,:-1].min()
        maxval = data[:,:-1].max()
        if subfig==1:
                cm = 'PuOr'
                p = ax1.pcolor(data, cmap=cm, alpha=0.8, edgecolors='w', linewidth=0.1, vmin=minval, vmax=maxval)
        else:
                cm = 'PuOr_r'
                p = ax1.pcolor(data, cmap=cm, norm=norm, alpha=0.8, edgecolors='w', linewidth=0.1, vmin=minval, vmax=maxval)
        # put the major ticks at the middle of each cell
        ax1.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
        ax1.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

        # want a more natural, table-like display
        ax1.invert_yaxis()
        ax1.xaxis.tick_top()
        ax1.set_frame_on(False)

        #ax1.set_title(title)
        if subfig == 1: ax1.set_xticklabels(labelsx, minor=False, rotation=90, fontsize=8)
        else: ax1.set_xticklabels([], minor=False)
        ax1.set_yticklabels(labelsy, minor=False, fontsize=8)


        # Remove ticks everywhere
        for t in ax1.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in ax1.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False

        # Add value to Total Score column
        for y in range(data.shape[0]):
                x = data.shape[1]-1
                plt.text(x + 0.5, y + 0.5, '%.1f' % data[y, x],
                horizontalalignment='center',
                verticalalignment='center', fontsize=8)

        cbar=plt.colorbar(p)
        cbar.ax.tick_params(labelsize=8)

# Add total score label to types
#phtypes = [p.replace('&','&&') for p in phtypes]
phtypes += ['TotalScore']
distlabels = np.array(results_distance.keys())[dorder].tolist()
enelabels = np.array(results_energy.keys())[eorder].tolist()
addheatmap(fig, 2, 1, dist_data, phtypes, distlabels, 'Distance scoring')
addheatmap(fig, 2, 2, ene_data, phtypes, enelabels, 'Energy instide ph4')
#ax.set_title('Scoring MDMix probes by pharmacophore overlapping')
plt.tight_layout()
plt.savefig(outname)

