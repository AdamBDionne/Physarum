using LightGraphs 

include("fxns/compute_phase_diagrams.jl")
include("fxns/phase_diagram_helper_functions.jl")
include("fxns/model.jl")
include("fxns/model_helper_functions.jl")

#Physarum 
    g, avg_radii, Len = import_MATLAB_network("sample_spring_two", false)
    params = graph_precomputation(g, avg_radii, Len)

#Example graphs
    # n = 21
    # g = LightGraphs.grid([n, 1]; periodic=false)
    # params = graph_precomputation(g, [], [])

# Save parameters
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b, modes = params

                                    #ϵ_s     #ϵ_c     #size
@time result = constructPhaseDiagram(g, 0.28, 0.45, 0.0, 1.5, 100)
