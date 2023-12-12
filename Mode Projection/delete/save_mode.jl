using DelimitedFiles 
mode_num = Ne-1

destination = "Raw Data/sample_spring_three/sample_spring_three"

g, avg_radii, Len = import_MATLAB_network("sample_spring_three", true)
params = graph_precomputation(g, avg_radii, Len)

# Save parameters
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, R0, Cq, mu, b = params


modes = copy(elastic_modes)
for j = 1:Ne
    modes[:,j] = elastic_modes[:,j] .- mean(elastic_modes[:,j])
    modes[:,j] = modes[:,j] ./ maximum(abs.(modes[:,j]))  
end

output = string(destination,"_MODES.txt")
open(output,"w") do io
    writedlm(io,modes)
end


# output = string(destination,"_mode_",mode_num,".txt")
# open(output,"w") do io
#     writedlm(io,modes[:,mode_num])
# end

