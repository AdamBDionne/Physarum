# This script saves the data collected from a simulation
using DelimitedFiles 
include("fxns/model.jl")

destination = "Data/Sample one/Sample one"

# Save parameters
params = [α[1];β[1];kappa[1];ϵ_s[1];ϵ_c[1];D[1]]

output = string(destination,"_PARAMETERS.txt")
open(output,"w") do io
    writedlm(io,["params";params;"init";u0])
end


# Save graph edges
    # edgesM = zeros(Int,ne(g),2)
    # for (i,e) in enumerate(edges(g))
    #     edgesM[i,1] = src(e)
    #     edgesM[i,2] = dst(e)
    # end

    # output = string(destination, "_EDGES.txt")
    # open(output,"w") do io
    #     writedlm(io,edgesM)
    # end

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
    #     #n = 5
    #     nodes = zeros(Int,nv(g),3)
    #     for i in 1:nv(g)
    #         nodes[i,1] = i 
    #         nodes[i,2] = mod(i-1,n)
    #         nodes[i,3] = floor((i-1)/n)
    #     end

# output = string(destination, "_NODES.txt")
# open(output,"w") do io
#     writedlm(io,nodes)
# end



# Save solution for animation
numFrames = fps*(t_f-t_i)
odeSol = [sol(t_i+(i/fps)) for i in 1:numFrames]

temp = zeros(numFrames,2Ne)
for i = 1:numFrames
    temp[i,:] = odeSol[i]
end

#Find R/avg{R} 
for i = 1:Ne
    curr_sig = temp[:,i]
    curr_sig = sqrt.(curr_sig ./ (pi .* Len[i]))
    temp[:,i] = curr_sig  ./ mean(curr_sig)
end
#Find C/avg{C}
for i = Ne+1:2Ne
    curr_sig = temp[:,i]
    temp[:,i] = curr_sig  ./ mean(curr_sig)
end
for i = 1:numFrames
    odeSol[i] = temp[i,:]
end

output = string(destination, "_ODE_SOL.txt")
open(output,"w") do io
    writedlm(io,odeSol)
end

# Save flows
flows = [findFlows(sol(t_i+(i/fps))) for i in 1:fps*(t_f-t_i)]

output = string(destination, "_FLOWS.txt")
open(output,"w") do io
    writedlm(io,flows)
end

println("SAVED!")