import numpy as npy
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

a= npy.loadtxt("/tmp/mdmix_orient_hpR8sO")
result = PCA(a)
x = []
y = []
z = []
for item in result.Y:
 x.append(item[0])
 y.append(item[1])
 z.append(item[2])

plt.close('all') # close all latent plotting windows
fig1 = plt.figure() # Make a plotting figure
#ax = Axes3D(fig1) # use the plotting figure to create a Axis3D object.
pltData = [x,y,z] 
colors=npy.linspace(0,256,len(x))


plt.subplot(131)
plt.scatter(pltData[0], pltData[1], marker='o', c=colors) # make a scatter plot of blue dots from the data
 
# make simple, bare axis lines through space:
#xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis 
#ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
#yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
#ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
#zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
#ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.
 
# label the axes 
plt.xlabel("PC1: %.3f"%result.fracs[0]) 
plt.ylabel("PC2: %.3f"%result.fracs[1])
#ax.set_zlabel("z-axis label")
#plt.title("The title of the plot")
#plt.show() # show the plot

plt.subplot(132)
plt.scatter(pltData[0], pltData[2], marker='o', c=colors)
plt.xlabel("PC1: %.3f"%result.fracs[0])
plt.ylabel("PC3: %.3f"%result.fracs[2])

plt.subplot(133)
plt.scatter(pltData[1], pltData[2], marker='o', c=colors)
plt.xlabel("PC2: %.3f"%result.fracs[1])
plt.ylabel("PC3: %.3f"%result.fracs[2])

plt.savefig("pca_plot.png")
