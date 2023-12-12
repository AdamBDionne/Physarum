# Adam Dionne

# code to animate nutrient transport

#-------- IMPORT PACKAGES -------#
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as manimation
import math
from os.path import exists
import numpy

#-- TO SPECIALIZE
file_suffix = "Raw Data/grid_5x5/grid_5x5"

destination = file_suffix
source = file_suffix
title = "Nutrient Transport"

still_frame = True

#--

#IMPORT DATA
with open(source + '_EDGES.txt') as f:
    edgeList = [[int(x) for x in line.split()] for line in f]

with open(source + '_NODES.txt') as f:
    nodeList = [[float(x) for x in line.split()] for line in f]

with open(source + '_ODE_SOL.txt') as f:
    sol = [[float(x) for x in line.split()] for line in f]

with open(source + '_FLOWS_IN.txt') as f:
    arrows_in = [[float(x) for x in line.split()] for line in f]

with open(source + '_FLOWS_OUT.txt') as f:
    arrows_out = [[float(x) for x in line.split()] for line in f]

with open(source + '_LEN.txt') as f:
    edge_lens = [[float(x) for x in line.split()][0] for line in f]



#SET FIGURE SIZE
fig, ax = plt.subplots(figsize=(7, 3))

#initialize graph
g = nx.Graph()
init = sol[0]
ne = len(edgeList)
NumPart = int((len(sol[0])-2*ne)/2)

#find relevant rescalings for volume 
maxV = max([max(i[0:ne-1]) for i in sol])
minV = min([min(i[0:ne-1]) for i in sol])

avg = (maxV+minV)/2
dist = max([maxV-avg,avg-minV])

scale = 5/(dist)
shift = 4

#find scaling for color bar
maxC = max([max(i[ne:2*ne-1]) for i in sol])
minC = min([min(i[ne:2*ne-1]) for i in sol])

#rescale flows
maxQ = 0
for x in range(len(arrows_in)):
    currMax = abs(max(arrows_in[x], key=abs))
    if currMax > maxQ:
        maxQ = currMax
#maxQ = max(max(arrows, key=abs), key=abs)
minQ = min(min(arrows_in))
#maxxQ = max([maxQ,abs(minQ)])
rescale = abs(1/maxQ)

#add nodes
for i in nodeList:
    g.add_node(i[0],pos=(i[1],i[2]))


#arrays with middle positions of every edge
midX_IN = [0 for i in range(ne)]
midY_IN = [0 for i in range(ne)]
midX_OUT = [0 for i in range(ne)]
midY_OUT = [0 for i in range(ne)]

#vector from node (j) to node (i)
deltaX = [0 for i in range(ne)]
deltaY = [0 for i in range(ne)]

#counter
c = 0 

#add edges & find positions for flow arrows 

for i in edgeList:
    g.add_edge(i[0],i[1])
    #to place arrows, we find the mid point of every edge
    # to add: angles
    pos1 = g.nodes[i[0]]["pos"]
    pos2 = g.nodes[i[1]]["pos"]

    #enters at node (j: pos2)
    midX_IN[c] = 0.25*pos1[0] + 0.75*pos2[0]
    midY_IN[c] = 0.25*pos1[1] + 0.75*pos2[1]

    #exits at node (i: pos1)
    midX_OUT[c] = 0.75*pos1[0] + 0.25*pos2[0]
    midY_OUT[c] = 0.75*pos1[1] + 0.25*pos2[1]

    #point from node (j: pos2) to node (i: pos1)
    deltaX[c] = pos1[0]-pos2[0]
    deltaY[c] = pos1[1]-pos2[1]
    c += 1



#get positions
pos = nx.get_node_attributes(g,'pos')

#keep track of progress
perc = 1 
tscale = len(sol)/10

def update(i):
    #report progress
    global perc
    if perc*tscale < i:
        print(str(perc) + '0% DONE')
        perc += 1

    if i == len(sol):
        print("SAVING...")

    #clear plot and construct frame
    plt.clf()
    ax = plt.gca()
    ax.set_title(title)
    plt.axis('equal')
    plt.axis('off')

    #update concentrations and volumes
    curr = sol[i]
    colors = [curr[ne+i] for i in range(ne)]
    vols = [scale*(curr[i]-minV)+shift for i in range(ne)]


    nodes = nx.draw_networkx_nodes(g,pos,node_size=0)
    edges = nx.draw_networkx_edges(
        g,
        pos,
        node_size=0,
        edge_color= colors,
        edge_cmap= plt.cm.viridis,
        width=vols,
        edge_vmin = minC,
        edge_vmax = maxC,
        ax = ax
    )

    #color bar
    cbar = plt.colorbar(edges)
    cbar.ax.set_ylabel("Concentration")

    #add flow arrows (IN)
    Qx = [0.5*max(0.5,abs(rescale*arrows_in[i][j]))*numpy.sign(arrows_in[i][j])*deltaX[j] for j in range(ne)]
    Qy = [0.5*max(0.5,abs(rescale*arrows_in[i][j]))*numpy.sign(arrows_in[i][j])*deltaY[j] for j in range(ne)]

    plt.quiver(midX_IN, midY_IN, Qx, Qy,scale=1,units='xy', edgecolor='azure', facecolor='tomato',linewidth = 0.4, pivot="mid")

    Qx = [0.5*max(0.5,abs(rescale*arrows_out[i][j]))*numpy.sign(arrows_out[i][j])*deltaX[j] for j in range(ne)]
    Qy = [0.5*max(0.5,abs(rescale*arrows_out[i][j]))*numpy.sign(arrows_out[i][j])*deltaY[j] for j in range(ne)]

    plt.quiver(midX_OUT, midY_OUT, Qx, Qy,scale=1,units='xy', edgecolor='azure', facecolor='tomato',linewidth = 0.4, pivot="mid")

    #draw nutrient particles 
    for part_ind in range(NumPart):
        #find position
        edge_ind = round(sol[i][2*ne+NumPart+part_ind])
        #make sure we aren't using an interpolated indice
        if abs(edge_ind - sol[i][2*ne+NumPart+part_ind]) > 0.01:
            #use the next valid indice
            c = i
            decrease = False
            while abs(round(sol[c][2*ne+NumPart+part_ind]) - sol[c][2*ne+NumPart+part_ind]) > 0.01:
                if decrease:
                    c -= 1
                if c < len(sol)-1:
                    c+=1
                #handle edge case where at end of sol. need to go backwards instead
                else: 
                    decrease = True
                    c = i 
                    c -= 1
            edge_ind = round(sol[c][2*ne+NumPart+part_ind])
        edge_nodes = edgeList[edge_ind-1]
        node_i = pos[edge_nodes[0]]
        node_j = pos[edge_nodes[1]]
        edge_len = (((node_i[0]-node_j[0])**2)+((node_i[1]-node_j[1])**2))**(0.5)

        dist_along_vessel = sol[i][2*ne+part_ind] / edge_lens[edge_ind - 1]
        prop_along = dist_along_vessel/edge_len

        x_pos = node_j[0]*(1-prop_along) + node_i[0]*prop_along
        y_pos = node_j[1]*(1-prop_along) + node_i[1]*prop_along

        ax.scatter([x_pos], [y_pos], s=80, color="k", edgecolors="w", zorder=10000)

#

print("ANIMATING...")
fig = plt.gcf()
if still_frame:
    ani = FuncAnimation(fig, update, frames=1)
    saveas = destination +'_animation.gif'
else:
    ani = FuncAnimation(fig, update, frames=len(sol))
    saveas = destination +'_animation.mp4'

cont = True
if exists(saveas):
    cont = input("Did you want to overwrite this file?")

if cont:
    writervideo = manimation.FFMpegWriter(fps=60)
    #writervideo = manimation.PillowWriter(fps=60)
    #set DPI for quality. 100 is low res for testing, 1000 is to share 
    ani.save(saveas, writer=writervideo, dpi = 150)
    print('DONE!')