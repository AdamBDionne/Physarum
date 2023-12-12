#code to generate m,n hexagonal lattice

#-------- IMPORT PACKAGES -------#
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib as mpl

g = nx.hexagonal_lattice_graph(5, 5, periodic = False, with_positions = True)

nodeList = list(g.nodes())
edgeList = list(g.edges())
len(edgeList)
pos = nx.get_node_attributes(g, 'pos')

#save 
with open('hexagon_EDGES.txt', 'w') as f:
    for i in list(g.edges()):
        f.write(str( nodeList.index(i[0])+1 ))
        f.write('\t')
        f.write(str( nodeList.index(i[1])+1 ))
        f.write('\n')

count = 1
with open('hexagon_NODES.txt', 'w') as f:
    for i in nodeList:
        f.write(str(count))
        f.write('\t')
        f.write(str(pos[i][0]))
        f.write('\t')
        f.write(str(pos[i][1]))
        f.write('\n')
        count += 1