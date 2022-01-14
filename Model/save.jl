# Adam Dionne
# This script saves the data collected from a simulation
using DelimitedFiles 
include("fxns/model.jl")

destination = "grid_5x5"

# Sampling
t_i = 0
t_f = 60
fps = 60

# Save parameters
params = [α;β;kappa;ϵ_s;ϵ_c;D]

output0 = string(destination,"_PARAMETERS.txt")
open(output0,"w") do io
    writedlm(io,["params";params;"init";u0])
end


# Save graph edges
    edgesM = zeros(Int,ne(g),2)
    for (i,e) in enumerate(edges(g))
        edgesM[i,1] = src(e)
        edgesM[i,2] = dst(e)
    end

    output1 = string(destination, "_EDGES.txt")
    open(output1,"w") do io
        writedlm(io,edgesM)
    end

# Save nodes and their positions 
# Make sure to change depending on the graph being used 

    #PATH GRAPH

    # nodes = zeros(Int,nv(g),3)
    # for i in 1:nv(g)
    #     nodes[i,1] = i 
    #     nodes[i,2] = i 
    #     nodes[i,3] = 0
    # end

    # GRID GRAPH
    #n = 5
    nodes = zeros(Int,nv(g),3)
    for i in 1:nv(g)
        nodes[i,1] = i 
        nodes[i,2] = mod(i-1,n)
        nodes[i,3] = floor((i-1)/n)
    end

    output2 = string(destination, "_NODES.txt")
    open(output2,"w") do io
        writedlm(io,nodes)
    end


# Save solution for animation

odeSol = [sol(t_i+(i/fps)) for i in 1:fps*(t_f-t_i)]

output3 = string(destination, "_ODE_SOL.txt")
open(output3,"w") do io
    writedlm(io,odeSol)
end

# Save flows
flows = [findFlows(sol(t_i+(i/fps))) for i in 1:fps*(t_f-t_i)]

output4 = string(destination, "_FLOWS.txt")
open(output4,"w") do io
    writedlm(io,flows)
end

println("SAVED!")