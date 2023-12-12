# Animate model solution for path graphs / grid graphs 

#-------- IMPORT PACKAGES -------#
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as manimation
import numpy
from os.path import exists

#-- TO SPECIALIZE
file_suffix = "Path Graph"

destination = "Movies/" + file_suffix + "/" + file_suffix
source = "Data/" + file_suffix + "/" + file_suffix
title = "Path Graph Peristalsis"

still_frame = False
#--

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
fig, ax = plt.subplots(figsize=(8, 5))

#initialize graph
g = nx.Graph()
init = sol[0]
ne = len(edgeList)

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

    #title
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
    cbar.ax.set_ylabel("Concentration (AU)")

    #add flow arrows
    Qx = [min(1.0,3*abs(rescale*arrows[i][j]))*numpy.sign(arrows[i][j])*deltaX[j] for j in range(ne)]
    Qy = [min(1.0,3*abs(rescale*arrows[i][j]))*numpy.sign(arrows[i][j])*deltaY[j] for j in range(ne)]

    plt.quiver(midX, midY, Qx, Qy,scale=1,units='xy', edgecolor='azure', facecolor='tomato',linewidth = 0.4, pivot="mid")

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
    writervideo = manimation.FFMpegWriter(fps=30)
    #writervideo = manimation.PillowWriter(fps=60)
    #set DPI for quality, use low when testing 
    ani.save(saveas, writer=writervideo, dpi = 50)
    print('DONE!')