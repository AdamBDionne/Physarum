#requires a solution from simulate.jl 
#calculates the average displacement of particles over solution 
using SimpleWeightedGraphs
using PyPlot
using PyCall
using Polynomials 
using LightGraphs
using DelimitedFiles


function find_average_displacement(pos_time, ind_time)
    g_extra = SimpleWeightedGraph(Nv+2)
    for i = 1:Ne
        SimpleWeightedGraphs.add_edge!(g_extra,edgelist[i,1],edgelist[i,2],Len[i])
    end

    average_displacement = zeros(size(pos_time,1))

    #go over all time
    for j = 1:size(pos_time,2)
        curr_disp = zeros(size(pos_time,1),1)
        last_ind = 0 
        last_disp = 0
        #find displacement for each particle
        for k = 1:size(pos_time,1)
            start_pos = pos_0[j]
            pos = pos_time[k,j]

            start_ind = ind_0[j]
            edge_ind = ind_time[k,j]
            gt = false

            if edge_ind != last_ind
                #add nodes where the particle began and end
                traversal_graph = copy(g_extra) 
                SimpleWeightedGraphs.rem_edge!(traversal_graph, edgelist[start_ind,1], edgelist[start_ind,2])
                SimpleWeightedGraphs.rem_edge!(traversal_graph, edgelist[edge_ind,1], edgelist[edge_ind,2])

                SimpleWeightedGraphs.add_edge!(traversal_graph, edgelist[start_ind,1],Nv+1,Len[start_ind]-start_pos)
                SimpleWeightedGraphs.add_edge!(traversal_graph, Nv+1,edgelist[start_ind,2],start_pos)

                SimpleWeightedGraphs.add_edge!(traversal_graph, edgelist[edge_ind,1],Nv+2,Len[edge_ind]-pos)
                SimpleWeightedGraphs.add_edge!(traversal_graph, Nv+2,edgelist[edge_ind,2],pos)


                #find shortest path
                short_path = SimpleWeightedGraphs.enumerate_paths(SimpleWeightedGraphs.dijkstra_shortest_paths(traversal_graph, Nv+1), Nv+2)
                weights = SimpleWeightedGraphs.weights(traversal_graph)
                #find displacement
                disp = 0
                num_trav = length(short_path)
                for i = 1:num_trav-1
                    disp += weights[short_path[i],short_path[i+1]]
                end
                curr_disp[k] = disp
                last_ind = edge_ind
                if weights[short_path[num_trav-1],short_path[end]] == pos
                    last_disp = disp - pos
                    gt = false
                else
                    last_disp = disp - (Len[edge_ind] - pos)
                    gt = true
                end
                last_disp = disp - pos
            else
                if gt
                    curr_disp[k] = last_disp + (Len[edge_ind]-pos)
                else
                    curr_disp[k] = last_disp + pos
                end
                
            end
        
        end
        average_displacement = average_displacement .+ curr_disp
    end
    average_displacement = average_displacement ./ num_parts 
    return average_displacement
end

@time t = find_average_displacement(pos_time, ind_time)