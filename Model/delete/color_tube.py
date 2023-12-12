# Animate model solution

#-------- IMPORT PACKAGES -------#
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as manimation
import math
import numpy
from os.path import exists
import matplotlib as mpl
import numpy as np

#-- TO SPECIALIZE
file_suffix = "path_graph"

ind = 80

destination = "Processed Data/" + file_suffix + "/" + file_suffix
source = "Raw Data/" + file_suffix + "/" + file_suffix
title = "Nutrient Transport"



#IMPORT DATA
with open(source + '_EDGES.txt') as f:
    edgeList = [[int(x) for x in line.split()] for line in f]

with open(source + '_NODES.txt') as f:
    nodeList = [[float(x) for x in line.split()] for line in f]

with open(source + '_ODE_SOL.txt') as f:
    sol = [[float(x) for x in line.split()] for line in f]

with open(source + '_FLOWS.txt') as f:
    arrows = [[float(x) for x in line.split()] for line in f]

#SET FIGURE SIZE
fig, ax = plt.subplots(figsize=(7, 3))

#initialize graph
g = nx.Graph()
init = sol[0]
ne = len(edgeList)

#find relevant rescalings for volume 
maxV = max([max(i[0:ne-1]) for i in sol])
minV = min([min(i[0:ne-1]) for i in sol])

avg = (maxV+minV)/2
dist = max([maxV-avg,avg-minV])

scale = 20/(dist)
shift = 10

#find scaling for color bar
maxC = max([max(i[ne:2*ne-1]) for i in sol])
minC = min([min(i[ne:2*ne-1]) for i in sol])

#rescale flows
maxQ = 0
for x in range(len(arrows)):
    currMax = abs(max(arrows[x], key=abs))
    if currMax > maxQ:
        maxQ = currMax
#maxQ = max(max(arrows, key=abs), key=abs)
minQ = min(min(arrows))
#maxxQ = max([maxQ,abs(minQ)])
rescale = abs(1/maxQ)

#add nodes
for i in nodeList:
    g.add_node(i[0],pos=(i[1],i[2]))

#add edges & find positions for flow arrows 

#arrays with middle positions of every edge
midX = [0 for i in range(ne)]
midY = [0 for i in range(ne)]

#vector from node (j) to node (i)
deltaX = [0 for i in range(ne)]
deltaY = [0 for i in range(ne)]

#counter
c = 0 
for i in edgeList:
    g.add_edge(i[0],i[1])
    #to place arrows, we find the mid point of every edge
    # to add: angles
    pos1 = g.nodes[i[0]]["pos"]
    pos2 = g.nodes[i[1]]["pos"]
    midX[c] = 0.5*(pos1[0] + pos2[0])
    midY[c] = 0.5*(pos1[1] + pos2[1])
    #point from node (j: pos2) to node (i: pos1)
    deltaX[c] = pos1[0]-pos2[0]
    deltaY[c] = pos1[1]-pos2[1]
    c += 1

#get positions
pos = nx.get_node_attributes(g,'pos')

#keep track of progress
perc = 1 
tscale = len(sol)/10


i = ind

cmap = plt.cm.Spectral  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]

# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)

# define the bins and normalize
bounds = np.linspace(1, ne+1, ne+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#clear plot and construct frame
plt.clf()
ax = plt.gca()
ax.set_title(title)
plt.axis('equal')
plt.axis('off')

#update concentrations and volumes
curr = sol[i]
colors = [(i+1)/ne for i in range(ne)]
vols = [1+shift for i in range(ne)]

#title
nodes = nx.draw_networkx_nodes(g,pos,node_size=0)
edges = nx.draw_networkx_edges(
    g,
    pos,
    node_size=0,
    edge_color= colors,
    edge_cmap= cmap,
    width=vols,
    edge_vmin = 0,
    edge_vmax = 1,
    ax = ax 
)

#color bar
cbar = plt.colorbar(edges, location ='bottom', cmap=cmap, norm=norm,
    spacing='proportional', ticks = bounds, boundaries=bounds)
cbar.ax.set_ylabel(i)

#add flow arrows
Qx = [min(1.0,3*abs(rescale*arrows[i][j]))*numpy.sign(arrows[i][j])*deltaX[j] for j in range(ne)]
Qy = [min(1.0,3*abs(rescale*arrows[i][j]))*numpy.sign(arrows[i][j])*deltaY[j] for j in range(ne)]

#plt.quiver(midX, midY, Qx, Qy,scale=1,units='xy', edgecolor='azure', facecolor='tomato',linewidth = 0.4, pivot="mid")

saveas = destination + '_' + str(ind) + '_still.eps'
plt.savefig(saveas, format='eps')
print('DONE!')