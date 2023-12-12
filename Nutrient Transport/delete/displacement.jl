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
    for j  = 1:size(pos_time,1)
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

                #opt1

                #opt2 

                #opt3 

                #opt4

                #find shortest path
                short_path = SimpleWeightedGraphs.enumerate_paths(SimpleWeightedGraphs.dijkstra_shortest_paths(traversal_graph, Nv+1), Nv+2)
                weights = SimpleWeightedGraphs.weights(traversal_graph)
                #find displacement
                disp = 0
                for i = 1:length(short_path)-1
                    disp += weights[short_path[i],short_path[i+1]]
                end
                curr_avg += disp^2
            end
        end
        average_displacement[j] = sqrt.(curr_avg / num_parts)
    end
    return average_displacement
end

#@time r = find_average_displacement(pos_time, ind_time)

# #plot results 
# #scale to micrometers
# PyPlot.close()
# average_displacement_unitfull = average_displacement * Char_L * 10^3 
# rcParams = PyDict(matplotlib["rcParams"])
# rcParams["font.size"] = 18
# fig, ax = PyPlot.subplots(figsize=(13,9))
# plt.subplot(1,2,1)
# #time = [(10/60)*i for i in 1:size(pos_time,1)]
# time = [(1/60)*i*10/fps for i in 1:(1/2*fps)*(t_f-t_i)]
# PyPlot.plot(time, average_displacement_unitfull)
# PyPlot.xlabel("Time (minutes)")
# PyPlot.ylabel("Average Displacement <x> (cm)")

# plt.subplot(1,2,2)
# cut = 20
# PyPlot.loglog(time[cut:end], average_displacement_unitfull[cut:end])
# PyPlot.xlabel("Log[Time (minutes)]")
# PyPlot.ylabel("Log[Average Displacement <x> (cm)]")

# line_fit = Polynomials.fit(log.(time[cut:end]),log.(average_displacement_unitfull[cut:end]),1)
# PyPlot.plot(exp.(log.(time[cut:end])), exp.(line_fit.(log.(time[cut:end]))))
# PyPlot.title(string("Slope = ",line_fit[1]))


# plt.tight_layout()
# PyPlot.savefig("disp_modes.png")


# output = string("avg_disp_modes.txt")
# open(output,"w") do io
#     writedlm(io,average_displacement)
# end

