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

    j = size(pos_time,1)

    curr_avg = 0 
    for i = 1:size(pos_time,2)
        start_pos = pos_0[i]
        pos = pos_time[j,i]

        start_ind = ind_0[i]
        edge_ind = ind_time[j,i]
    
        if start_ind == edge_ind
            curr_avg += abs(start_pos - pos)
        else
                
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
            for i = 1:length(short_path)-1
                disp += weights[short_path[i],short_path[i+1]]
            end
            curr_avg += disp
        end
    end

    first_max = curr_avg / num_parts
    last_max = first_max - 1
    curr_max = first_max 
    counter = 1
    while curr_max > last_max
        last_max = curr_max
        j = size(pos_time,1) - counter
        curr_avg = 0 
        for i = 1:size(pos_time,2)
            start_pos = pos_0[i]
            pos = pos_time[j,i]
    
            start_ind = ind_0[i]
            edge_ind = ind_time[j,i]
        
            if start_ind == edge_ind
                curr_avg += abs(start_pos - pos)
            else
                    
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
                for i = 1:length(short_path)-1
                    disp += weights[short_path[i],short_path[i+1]]
                end
                curr_avg += disp
            end
        end
        curr_max = curr_avg / num_parts
        counter = counter + 1
    end
    return last_max
end
