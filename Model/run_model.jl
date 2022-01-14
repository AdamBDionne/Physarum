# Simulate model 

# Packages
using LightGraphs 
using DifferentialEquations
using Plots

# Helper code 
include("fxns/model.jl")
include("fxns/model_helper_functions.jl")

# Intialize graph (if needed)
n = 8
g = LightGraphs.grid([n, n]; periodic=false)

# Load parameters
        # For grids
         params = graph_precomputation(g, [], [])
        # For Physarum
        # params = graph_precomputation(g, [], [])

# Save parameters
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, ϵ_s, ϵ_c, surface_area, radius_sq = params

# Initial conditions - add some noise to concentrations 
#u0 = [V0;C0 + 0.01randn(Ne)]

# Solve model - to remove transient 
T = 100
tspan = (0.0, T)
prob = ODEProblem(h, u0, tspan)
@time sol = solve(prob, RK4())#, saveat = 950:0.01:1000)

#Solve model w/o transient
T = 60
tspan = (0.0, T)
prob = ODEProblem(h, sol(99), tspan)
@time sol = solve(prob, RK4())