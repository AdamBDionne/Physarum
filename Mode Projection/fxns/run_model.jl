include("model_helper_functions.jl")
function run_model(file_suffix, p, fps)
    global g; global avg_radii; global Len
    g, avg_radii, Len = import_MATLAB_network(file_suffix, true)


    params = graph_precomputation(g, avg_radii, Len)


    global Ne;  global Nv; global avg_radii; global Len; global V0; global kappa; global α; global D; global β; global Deff; global Dt; global τ; global C0; global E; global B; global K; global k; global L; global Ldag; global Mq_inv; global Minv1; global sum_Minv1; global elastic_modes; global elastic_eigvals; global edgelist; global out_neighbors; global in_neighbors; global Ein; global Eout; global EinT; global EoutT; global surface_area; global radius_sq; global Char_L; global Char_t; global R0; global Cq; global mu; global b 
    Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b = params

    global ϵ_s; global ϵ_c
    ϵ_s, ϵ_c = p

    shift = elastic_modes[:,Ne-1] + elastic_modes[:,Ne-2]
    u0 = [V0; V0 .+ (0.05 .*(shift ./ maximum(abs.(shift))))]
    #u0 = [V0 ; V0 .+ (0.05 .* rand(Ne))]
    #For projecting modes
    t_i = 150
    t_f = 150+480

    # t_i = 100
    # t_f = 150

    T = t_f
    tspan = (0.0, T)
    prob = ODEProblem(h, u0, tspan)
    sol = solve(prob, RK4(), saveat = t_i:1/fps:t_f)
    return sol, elastic_modes, Ne, u0
end