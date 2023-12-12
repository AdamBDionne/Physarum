# Simulate network model 

# Packages
using LightGraphs 
using DifferentialEquations
using Plots

# Helper code 
include("fxns/model.jl")
include("fxns/model_helper_functions.jl")

# Intialize graph 

        # Grid/path graph 
                # n = 13
                # g = LightGraphs.grid([n, 1]; periodic=false)
                # params = graph_precomputation(g, [], [])

        # Experimental Physarum data 
                file_suffix = "Sample one"
                g, avg_radii, Len = import_MATLAB_network(file_suffix, true)
                params = graph_precomputation(g, avg_radii, Len)

# Save parameters
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b, L = params

# Set ϵ_s, ϵ_c
p = [0.38 0.4]
ϵ_s, ϵ_c = p

# Initial Conditions

# Equilibrium + noise 
        #V0 = ones(Ne)
        #u0 = [V0;C0 + 0.01randn(Ne)]

# Two dominant modes phase shifted 

        shift = elastic_modes[:,Ne-1] #- elastic_modes[:,Ne-2]
        u0 = [V0; V0 .+ (0.05 .*(shift ./ maximum(abs.(shift))))]

#Solve model w/o transient
t_i = 150
t_f = 150+30
fps = 30
# t_i = 50
# t_f = 250
# fps = 100

tspan = (0.0, t_f)
prob = ODEProblem(h, u0, tspan)
@time sol = solve(prob, RK4(), saveat = t_i:1/fps:t_f)