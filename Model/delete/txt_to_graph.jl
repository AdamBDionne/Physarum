#Adam Dionne 

#This script takes .txt files from in our standard format
#and converts to a light graph g

using DelimitedFiles
using LightGraphs
using SimpleWeightedGraphs

filename = "Raw Data/hexagon_graph/hexagon_"

#-- read in edge and node file ---
nodes = readdlm(string(filename,"NODES.txt"))
g = SimpleGraph(size(nodes)[1])
edgeTXT = readdlm(string(filename,"EDGES.txt"))
for i in 1:size(edgeTXT)[1]
    add_edge!(g,Int(edgeTXT[i,1]),Int(edgeTXT[i,2]))
end


# #skipto might be useful in future

# data = readdlm("physarum.txt",'|')

# g = SimpleWeightedGraph()
# nodePos = zeros(Int,size(data)[1],2) 

# #make a dictionary of tuple -> node value
# #each time check if our node is already in it
# #if not, add into dictionary
# #use dict to add edge 
# function constructGraph()
#     g = SimpleWeightedGraph()
#     dict = Dict()
#     C = 1
#     currNode = 1
#     for i in 1:size(data)[1]
#         s = strip(data[i,1], ['(',')'])
#         s = split(s,",")
#         a = parse(Int,s[1])
#         b = parse(Int,s[2])
#         if haskey(dict,(a,b)) == false
#             dict[(a,b)] = C
#             C += 1
#             add_vertex!(g)
#         end
#         if data[i,2][1] != '{'
#             currNode = dict[(a,b)]
#         else
#             j = split(data[i,2],' ')
#             j = strip(j[4],[','])
#        #add_edge!(g,currNode,dict[(a,b)],parse(Float16,j))
#             add_edge!(g,currNode,dict[(a,b)],5)
#         end
#     end

#     h = SimpleWeightedGraph(nv(g))
#     for (i,e) in enumerate(edges(g))
#         a = inneighbors(g,src(e))
#         a = size(a)[1]
#         b = inneighbors(g,dst(e))
#         b = size(b)[1]
#         if a == 1 || b == 1
#             rem_edge!(g,dst(e),src(e))
#         end
#     end
#     for (i,e) in enumerate(edges(g))
#         a = inneighbors(g,src(e))
#         a = size(a)[1]
#         b = inneighbors(g,dst(e))
#         b = size(b)[1]
#         if a == 1 || b == 1
#             rem_edge!(g,dst(e),src(e))
#         end
#     end
#     for (i,e) in enumerate(edges(g))
#         if weight(e) != 0 
#             print(ne(h))
#             add_edge!(h,src(e),dst(e),weight(e))
#         end
#     end
    
#     return h, dict
# end



# #save graph in the format used for animateNetwork.py
# function saveGraph(filename,dict)
#     edgesM = zeros(Int,ne(g),2)
#     local cnt = 1
#     for (i,e) in enumerate(edges(g))
#         edgesM[cnt,1] = src(e)
#         edgesM[cnt,2] = dst(e)
#         cnt += 1
#     end

#     output1 = string(filename, "EDGES.txt")
#     open(output1,"w") do io
#         writedlm(io,edgesM)
#     end

#     nodesM = zeros(Int,length(dict),3)
#     cnt = 1
#     for (key, value) in dict
#         nodesM[cnt,1] = value
#         nodesM[cnt,2] = key[1]
#         nodesM[cnt,3] = key[2]
#         cnt += 1 
#     end

#     output2 = string(filename, "NODES.txt")
#     open(output2,"w") do io
#         writedlm(io,nodesM)
#     end
# end

# g, dict = constructGraph()
# saveGraph("physarum", dict)



    